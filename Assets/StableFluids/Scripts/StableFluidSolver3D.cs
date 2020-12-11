using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class StableFluidSolver3D : MonoBehaviour
{
    enum BCType
    {
        vertical,
        horizontal,
        depth,
        pressure
    }

    //SIMULATION PARAMETERS
    [SerializeField]
    bool StartSimulation = false;
    //renderer
    [SerializeField]
    FluidRenderer renderer;

    GameObject debugVelocity;

    //SOLVER PARAMETERS
    [SerializeField]
    int N = 128;
    [SerializeField]
    float dt = 0.2f;
    float timer = 0;

    //FLUID PARAMETERS
    [SerializeField]
    float dissipationRate = 0.2f;
    [SerializeField]
    float diffusionRate = 0f;
    [SerializeField]
    float viscosity = 0.000001f;

    //2D arrays for veloctiy and dens
    float[,,] dens;
    float[,,] densPrev;
    float[,,] velX;
    float[,,] velXPrev;
    float[,,] velY;
    float[,,] velYPrev;
    float[,,] velZ;
    float[,,] velZPrev;

    //handling mouse input
    Vector2 lastMousePoint;
    Vector2 mouseSpeed;
    [SerializeField]
    float mouseSpeedFactor = 1;

    private void Start()
    {
        //2D arrays for veloctiy and dens
        dens = new float[N, N, N];
        densPrev = new float[N, N, N];
        velX = new float[N, N, N];
        velXPrev = new float[N, N, N];
        velY = new float[N, N, N];
        velYPrev = new float[N, N, N];
        velZ = new float[N, N, N];
        velZPrev = new float[N, N, N];

        renderer.renderScalarField3D(ref dens);

        debugVelocity = new GameObject();
        debugVelocity.name = "debugVelocity";
    }

    //MAIN LOOP

    #region mainLoop

    void Vstep()
    {
        diffuse(velXPrev, velX, viscosity, dt); setBoundaryConditions(BCType.vertical, velX);
        diffuse(velYPrev, velY, viscosity, dt); setBoundaryConditions(BCType.horizontal, velY);
        diffuse(velZPrev, velZ, viscosity, dt); setBoundaryConditions(BCType.horizontal, velZ);
        //projection step
        project(velXPrev, velYPrev,velZPrev, velX, velY);
        //Advection step
        transport(velX, velXPrev, velXPrev, velYPrev, velZPrev); setBoundaryConditions(BCType.vertical, velX);
        transport(velY, velYPrev, velXPrev, velYPrev, velZPrev); setBoundaryConditions(BCType.horizontal, velY);
        transport(velZ, velZPrev, velXPrev, velYPrev, velZPrev); setBoundaryConditions(BCType.horizontal, velZ);
        //projection step -> solving Poisson equation
        project(velX, velY,velZPrev, velXPrev, velYPrev);
    }

    void Sstep()
    {
        //diffusion step
        diffuse(densPrev, dens, diffusionRate, dt);
        //advection step
        transport(dens, densPrev, velX, velY, velZ); setBoundaryConditions(BCType.pressure, dens);
        //dissipation step
        dissipate(dens);
    }

    private void Update()
    {
        if (StartSimulation)
        {
            timer += Time.deltaTime;
            if (timer > dt)
            {
                Vstep();
                Sstep();
                timer -= dt;
                renderer.renderScalarField3D(ref dens);
            }
        }
        mouseSpeed = new Vector2(Input.mousePosition.x, Input.mousePosition.y) - lastMousePoint;
        lastMousePoint = Input.mousePosition;
    }

    #endregion

    //SOLVER METHODS

    #region solverMethods

    void diffuse(float[,,] x, float[,,] x0, float diff, float dt)
    {
        float diffusionRate = dt * diff * N * N;
        NumericalMethods.GaussSeidelRelaxLinearSolver3D(x, x0, diffusionRate, 1 + 6 * diffusionRate, N);
    }

    void project(float[,,] Ux, float[,,] Uy, float[,,] Uz, float[,,] pressure, float[,,] nDiv)
    {
        //compute divergence in the grid
        for (int k = 1; k< N - 1; k++)
        {
            for (int j = 1; j < N - 1; j++)
            {
                for (int i = 1; i < N - 1; i++)
                {
                    //computing negative divergence which will be the right hand of the equation for pressure linear system
                    nDiv[i, j,k] = -NumericalMethods.discreteVelocityDivergence3D(Ux, Uy, Uz, i, j,k, N);
                    pressure[i, j,k] = 0;
                }
            }

        }
        //reset bondary conditions since we made changes in the grid
        setBoundaryConditions(BCType.pressure, nDiv); setBoundaryConditions(BCType.pressure, pressure);
        //solve Poisson equation using GaussianSeidel solver
        NumericalMethods.GaussSeidelRelaxLinearSolver3D(pressure, nDiv, 1, 6, N);
        //Hodge decomposition : any field = mass conserving + gradient
        //we want ONLY the mass conserving part (= initial field - gradient)
        for (int k = 1; k < N - 1; k++)
        {
            for (int j = 1; j < N - 1; j++)
            {
                for (int i = 1; i < N - 1; i++)
                {
                    Ux[i, j, k] -= NumericalMethods.centralDifferences3D(pressure, i, j,k, 1, (1 / (float)N));
                    Uy[i, j,k] -= NumericalMethods.centralDifferences3D(pressure, i, j,k, 0, (1 / (float)N));
                    Uz[i, j,k] -= NumericalMethods.centralDifferences3D(pressure, i, j,k, 2, (1 / (float)N));
                }
            }
        }
        setBoundaryConditions(BCType.vertical, Ux);
        setBoundaryConditions(BCType.horizontal, Uy);
        setBoundaryConditions(BCType.horizontal, Uz);
    }

    void transport(float[,,] s, float[,,] s0, float[,,] Ux, float[,,] Uy, float[,,] Uz)
    {
        for (int k = 1; k < N - 1; k++)
        {
            for (int i = 1; i < N - 1; i++)
            {
                for (int j = 1; j < N - 1; j++)
                {
                    Vector3 prevParticlePos = NumericalMethods.TraceParticle3D(i, j,k, Ux, Uy, Uz,dt,N);
                    //get the 6 neighbours cells around position (x,y)
                    //set back the interpolated value based on these 6 cells back at the initial grid cell (i,j,k)      
                    s[i, j, k] = NumericalMethods.trilinearInterpolation(s0, prevParticlePos, N);
                     
                }
            }
        }
    }


    void setBoundaryConditions(BCType bctype, float[,,] x)
    {
        int iBCfactor = 1;
        int jBCfactor = 1;
        int kBCfactor = 1;
        switch (bctype)
        {
            case BCType.vertical:
                iBCfactor = -1;
                jBCfactor = -1;
                kBCfactor = 1;
                break;
            case BCType.horizontal:
                iBCfactor = 1;
                jBCfactor = -1;
                kBCfactor = 1;
                break;
            case BCType.depth:
                iBCfactor = 1;
                jBCfactor = 1;
                kBCfactor = -1;
                break;
            case BCType.pressure:
                break;
        }
        //borders along the grid      
        for (int i = 1; i < N - 1; i++)
        {
            for (int k = 1; k < N - 1; k++)
            {
                x[i, 0, k] = jBCfactor * x[i, 1, k];
                x[i, N - 1, k] = jBCfactor*x[i, N - 2, k];
            }
        }        
        for (int j = 1; j < N - 1; j++)
        {
            for (int i = 1; i < N - 1; i++)
            {
                x[i, j, 0] = kBCfactor * x[i, j, 1];
                x[i, j, N - 1] = kBCfactor * x[i, j, N - 2];
            }
        }
        for (int k = 1; k< N - 1; k++)
        {
            for (int j = 0; j < N - 1; j++)
            {
                x[0, j, k] = iBCfactor * x[1, j, k];
                x[N - 1, j, k] = iBCfactor * x[N - 2, j, k];
            }
        }
        //special case for corners
        //BOTTOM FACE
        x[0, 0, 0] = (x[1, 0, 0] + x[0, 1, 0] + x[0, 0, 1]) / 3f; //BL corner
        x[0, N - 1, 0] = (x[1, N - 1, 0] + x[0, N - 2, 0] + x[0, N - 1, 1]) / 3f; //BR corner
        x[0, N - 1, N-1] = (x[1, N - 1, N-1] + x[0, N - 2, N-1] + x[0, N - 1, N-2]) / 3f; //TR corner
        x[0, 0, N-1] = (x[1, 0, N-1] + x[0, 1, N-1] + x[0, 0, N-2]) / 3f; //TL corner
        //TOP FACE
        x[N - 1, 0, 0] = (x[N - 2, 0, 0] + x[N - 1, 1, 0] + x[N - 1, 0, 1]) / 3f; //BL corner
        x[N - 1, N - 1, 0] = (x[N - 2, N - 1, 0] + x[N - 1, N - 2, 0] + x[N - 1, N - 1, 1]) / 3f; //BR corner
        x[N - 1, N - 1, N-1] = (x[N - 2, N - 1, N-1] + x[N - 1, N - 2, N-1] + x[N - 1, N - 1, N-2]) / 3f; //TR corner
        x[N - 1, 0, N-1] = (x[N - 2, 0, N-1] + x[N - 1, 1, N-1] + x[N - 1, 0, N-2]) / 3f; //TL corner

    }

    void dissipate(float[,,] s)
    {
        for (int k = 1; k < N; k++)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    s[i, j,k] *= (1 / (1 + dt * dissipationRate));
                }
            }
        }
    }
    #endregion    

    // USER CONTROL METHODS

    #region userControls

    public void addDensitySource()
    {
        Vector3Int mousePos = renderer.getMousePosOnGrid3D();
        addDensity(mousePos.x, mousePos.y,mousePos.z, 100);
        print("adding d");
    }

    public void addVelocitySource()
    {
        Vector3Int mousePos = renderer.getMousePosOnGrid3D();
        if (mousePos.x > 0 && mousePos.x < N && mousePos.y > 0 && mousePos.y < N)
        {
            addDensity(mousePos.x, mousePos.y,mousePos.z, 10);
            addVelocity(mousePos.x, mousePos.y,mousePos.z, mouseSpeed.x*mouseSpeedFactor, mouseSpeed.y * mouseSpeedFactor, 0);
        }
    }
    void addDensity(int i, int j,int k, float amount)
    {
        dens[i, j,k] += amount;
    }

    void addVelocity(int i, int j, int k, float amountX, float amountY, float amountZ)
    {
        this.velX[i, j, k] += amountX;
        this.velY[i, j, k] += amountY;
        this.velZ[i, j, k] += amountZ;
    }
    #endregion
}
