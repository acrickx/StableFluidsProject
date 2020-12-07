using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class StableFluidSolver3D : MonoBehaviour
{
    enum BCType
    {
        type1,
        type2,
        type3,
        type0
    }

    //SIMULATION PARAMETERS
    [SerializeField]
    bool StartSimulation = false;
    //renderer
    [SerializeField]
    SimpleRenderer renderer;

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
        diffuse(velXPrev, velX, viscosity, dt); setBoundaryConditions(BCType.type1, velX);
        diffuse(velYPrev, velY, viscosity, dt); setBoundaryConditions(BCType.type2, velY);
        diffuse(velZPrev, velZ, viscosity, dt); setBoundaryConditions(BCType.type2, velZ);
        //projection step
        project(velXPrev, velYPrev,velZPrev, velX, velY);
        //Advection step
        transport(velX, velXPrev, velXPrev, velYPrev, velZPrev); setBoundaryConditions(BCType.type1, velX);
        transport(velY, velYPrev, velXPrev, velYPrev, velZPrev); setBoundaryConditions(BCType.type2, velY);
        transport(velZ, velZPrev, velXPrev, velYPrev, velZPrev); setBoundaryConditions(BCType.type2, velZ);
        //projection step -> solving Poisson equation
        project(velX, velY,velZPrev, velXPrev, velYPrev);
    }

    void Sstep()
    {
        //diffusion step
        diffuse(densPrev, dens, diffusionRate, dt);
        //advection step
        transport(dens, densPrev, velX, velY, velZ); setBoundaryConditions(BCType.type0, dens);
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
        if (Input.GetKeyDown(KeyCode.D))
        {
            debugDens();
        }
    }

    #endregion

    //SOLVER METHODS

    #region solverMethods

    void diffuse(float[,,] x, float[,,] x0, float diff, float dt)
    {
        float a = dt * diff * (N-2) * (N-2);
        GaussSeidelRelaxLinearSolver(x, x0, a, 1 + 6 * a);
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
                    nDiv[i, j,k] = -discreteVelocityDivergence(Ux, Uy, Uz, i, j,k);
                    pressure[i, j,k] = 0;
                }
            }

        }
        //reset bondary conditions since we made changes in the grid
        setBoundaryConditions(BCType.type0, nDiv); setBoundaryConditions(BCType.type0, pressure);
        //solve Poisson equation using GaussianSeidel solver
        GaussSeidelRelaxLinearSolver(pressure, nDiv, 1, 6);
        //Hodge decomposition : any field = mass conserving + gradient
        //we want ONLY the mass conserving part (= initial field - gradient)
        for (int k = 1; k < N - 1; k++)
        {
            for (int j = 1; j < N - 1; j++)
            {
                for (int i = 1; i < N - 1; i++)
                {
                    Ux[i, j,k] -= centralDifferences(pressure, i, j,k, 0, (1 / (float)N));
                    Uy[i, j,k] -= centralDifferences(pressure, i, j,k, 1, (1 / (float)N));
                    Uz[i, j,k] -= centralDifferences(pressure, i, j,k, 2, (1 / (float)N));
                }
            }
        }
        setBoundaryConditions(BCType.type1, Ux);
        setBoundaryConditions(BCType.type2, Uy);
        setBoundaryConditions(BCType.type2, Uz);
    }

    void transport(float[,,] s, float[,,] s0, float[,,] Ux, float[,,] Uy, float[,,] Uz)
    {
        for (int k = 1; k < N - 1; k++)
        {
            for (int i = 1; i < N - 1; i++)
            {
                for (int j = 1; j < N - 1; j++)
                {
                    Vector3 prevParticlePos = TraceParticle(i, j,k, Ux, Uy, Uz);
                    //get the 4 neighbours cells around position (x,y)
                    int x0 = Mathf.Clamp(Mathf.FloorToInt(prevParticlePos.x), 1, N - 2);
                    int x1 = Mathf.Clamp(x0 + 1, 0, N - 1);
                    int y0 = Mathf.Clamp(Mathf.FloorToInt(prevParticlePos.y), 1, N - 2);
                    int y1 = Mathf.Clamp(y0 + 1, 0, N - 1);
                    int z0 = Mathf.Clamp(Mathf.FloorToInt(prevParticlePos.z), 1, N - 2);
                    int z1 = Mathf.Clamp(z0 + 1, 0, N - 1);
                    //get the deltas between x/x0 and y/y0 to weight the linear combination for the density
                    float wx1 = prevParticlePos.x - x0;
                    float wx0 = 1 - wx1;
                    float wy1 = prevParticlePos.y - y0;
                    float wy0 = 1 - wy1;
                    float wz1 = prevParticlePos.z - z0;
                    float wz0 = 1 - wz1;
                    //set back the interpolated value based on these 4 cells back at the initial grid cell (i,j)                
                    s[i, j, k] = wx0 * (wy0 * (wz0 * s0[x0, y0, z0] + wz1 * (s0[x0, y0, z1])) + wy1 * (wz0 * s0[x0, y1, z0] + wz1 * (s0[x0, y1, z1])))
                               + wx1 * (wy0 * (wz0 * s0[x1, y0, z0] + wz1 * (s0[x1, y0, z1])) + wy1 * (wz0 * s0[x1, y1, z0] + wz1 * (s0[x1, y1, z1])));                              
                }
            }
        }
    }

    //trace particle back to previous position
    Vector3 TraceParticle(int i, int j, int k, float[,,] Ux, float[,,] Uy, float[,,] Uz)
    {
        //tracing the particle backwards in time 
        float a = (i - (dt * N) * Ux[i, j,k]); //new vertical position from center of cell
        float b = (j - (dt * N) * Uy[i, j,k]); //new horizontal position
        float c = (k - (dt * N) * Uz[i, j,k]); //new depth position
        // make sure the position of the particle is inside the grid; if not set it at the limit
        a = Mathf.Clamp(a, 0.99f, N - 0.01f);
        b = Mathf.Clamp(b, 0.99f, N - 0.01f);
        c = Mathf.Clamp(c, 0.99f, N - 0.01f);
        // return traced position
        return new Vector3(a, b, c);
    }

    void setBoundaryConditions(BCType bctype, float[,,] x)
    {
        int iBCfactor = 1;
        int jBCfactor = 1;
        int kBCfactor = 1;
        switch (bctype)
        {
            case BCType.type1:
                iBCfactor = -1;
                jBCfactor = -1;
                kBCfactor = 1;
                break;
            case BCType.type2:
                iBCfactor = 1;
                jBCfactor = -1;
                kBCfactor = 1;
                break;
            case BCType.type3:
                iBCfactor = 1;
                jBCfactor = 1;
                kBCfactor = -1;
                break;
            case BCType.type0:
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
        //TOP FACE
        x[0, 0, 0] = (x[1, 0, 0] + x[0, 1, 0] + x[0, 0, 1]) / 3f; //BL corner
        x[0, N - 1, 0] = (x[1, N - 1, 0] + x[0, N - 2, 0] + x[0, N - 1, 1]) / 3f; //BR corner
        x[0, N - 1, N-1] = (x[1, N - 1, N-1] + x[0, N - 2, N-1] + x[0, N - 1, N-2]) / 3f; //TR corner
        x[0, 0, N-1] = (x[1, 0, N-1] + x[0, 1, N-1] + x[0, 0, N-2]) / 3f; //TL corner
        //BOTTOM FACE
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

    //NUMERICAL METHODS

    #region numericalMethods

    float centralDifferences(float[,,] x, int i, int j, int k, int direction, float dx)
    {
        if (direction == 0)
            return (x[i + 1, j,k] - x[i - 1, j,k]) / (2f * dx);
        else if(direction ==1)
            return (x[i, j + 1,k] - x[i, j - 1,k]) / (2f * dx);
        else
            return (x[i, j, k+1] - x[i, j - 1, k-1]) / (2f * dx);
    }

    void GaussSeidelRelaxLinearSolver(float[,,] x, float[,,] x0, float leftHand, float rightHand, int limitStep = 16)
    {
        for (int step = 0; step < limitStep; step++)
        {
            for (int k = 1; k < N - 1; k++)
            {
                for (int i = 1; i < N - 1; i++)
                {
                    for (int j = 1; j < N - 1; j++)
                    {
                        //linear combination of "cross" neighbours
                        x[i, j, k] = (x0[i, j, k] + leftHand * (x[i + 1, j, k]
                                                       + x[i - 1, j, k]
                                                       + x[i, j + 1, k]
                                                       + x[i, j - 1, k]
                                                       + x[i, j, k + 1]
                                                       + x[i, j, k - 1])) / (float)rightHand;
                    }
                }
            }
        }
    }

    float discreteVelocityDivergence(float[,,] Ux, float[,,] Uy, float[,,] Uz, int i, int j, int k)
    {
        return (centralDifferences(Ux, i, j,k, 0, N) + centralDifferences(Uy, i, j,k, 1, N)) + centralDifferences(Uz,i,j,k,2,N);
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

    void debugDens()
    {
        float[,] d = new float[N, N];
        for (int i = 1; i < N - 1; i++)
        {
            for (int j = 1; j < N - 1; j++)
            {
                d[i, j] = dens[i, j, 10];
            }
        }
        renderer.renderScalarField(ref d);
    }
    #endregion
}
