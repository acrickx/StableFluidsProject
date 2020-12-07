using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class StableFluidSolver : MonoBehaviour
{
    enum BCType
    {
        type1,
        type2,
        type0
    }

    //SIMULATION PARAMETERS
    [SerializeField]
    bool StartSimulation = false;
    //renderer
    [SerializeField]
    SimpleRenderer renderer;

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
    float[,] dens;
    float[,] densPrev;
    float[,] velX;
    float[,] velXPrev;
    float[,] velY;
    float[,] velYPrev;

    //handling mouse input
    Vector2 lastMousePoint;
    Vector2 mouseSpeed;
    [SerializeField]
    float mouseSpeedFactor = 1;

    //debugging
    GameObject debugVelocity;

    private void Start()
    {
        //2D arrays for veloctiy and dens
        dens = new float[N, N];
        densPrev = new float[N, N];
        velX = new float[N, N];
        velXPrev = new float[N, N];
        velY = new float[N, N];
        velYPrev = new float[N, N];

        debugVelocity = new GameObject();
        debugVelocity.name = "debugVelocity";
    }

    //MAIN LOOP
    #region mainLoop

    void Vstep()
    {
        diffuse(velXPrev, velX, viscosity); 
        setBoundaryConditions(BCType.type1, velX);
        diffuse(velYPrev, velY, viscosity); 
        setBoundaryConditions(BCType.type2, velY);
        //projection step
        project(velXPrev, velYPrev, velX, velY);
        //Advection step
        transport(velX, velXPrev, velXPrev, velYPrev); 
        setBoundaryConditions(BCType.type1, velX);
        transport(velY, velYPrev, velXPrev, velYPrev); 
        setBoundaryConditions(BCType.type2, velY);
        //projection step -> solving Poisson equation
        project(velX, velY, velXPrev, velYPrev);
    }

    void Sstep()
    {
        //diffusion step
        diffuse(densPrev, dens, diffusionRate);
        //advection step
        transport(dens, densPrev, velX, velY); 
        setBoundaryConditions(BCType.type0, dens);
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
                renderer.renderScalarField(ref dens);
            }
        }
        mouseSpeed = new Vector2(Input.mousePosition.x, Input.mousePosition.y) - lastMousePoint;
        lastMousePoint = Input.mousePosition;

        if (Input.GetKeyDown(KeyCode.V))
        {
            debugVelocities();
        }

        if (Input.GetKeyDown(KeyCode.Escape))
        {
            clearVelocities();
        }
    }

    #endregion

    //SOLVER METHODS

    #region solverMethods

    void diffuse(float[,] x, float[,] x0, float diff)
    {
        float a = dt * diff * N * N;
        GaussSeidelRelaxLinearSolver(x, x0, a, 1 + 4 * a);
    }

    void project(float[,] Ux, float[,] Uy, float[,] pressure, float[,] nDiv)
    {
        //compute divergence in the grid
        for (int j = 1; j < N - 1; j++)
        {
            for (int i = 1; i < N - 1; i++)
            {
                //computing negative divergence which will be the right hand of the equation for pressure linear system
                nDiv[i, j] = -discreteVelocityDivergence(Ux, Uy, i, j);
                pressure[i, j] = 0;
            }
        }
        //reset bondary conditions since we made changes in the grid
        setBoundaryConditions(BCType.type0, nDiv); setBoundaryConditions(BCType.type0, pressure);
        //solve Poisson equation using GaussianSeidel solver
        GaussSeidelRelaxLinearSolver(pressure, nDiv, 1, 4);
        //Hodge decomposition : any field = mass conserving + gradient
        //we want ONLY the mass conserving part (= initial field - gradient)
        for (int j = 1; j < N - 1; j++)
        {
            for (int i = 1; i < N - 1; i++)
            {
                Ux[i, j] -= centralDifferences(pressure, i, j, 0, (1 / (float)N));
                Uy[i, j] -= centralDifferences(pressure, i, j, 1, (1 / (float)N));
            }
        }
        setBoundaryConditions(BCType.type1, Ux);
        setBoundaryConditions(BCType.type2, Uy);
    }

    void transport(float[,] s, float[,] s0, float[,] Ux, float[,] Uy)
    {
        for (int i = 1; i < N - 1; i++)
        {
            for (int j = 1; j < N - 1; j++)
            {
                Vector2 prevParticlePos = TraceParticle(i, j, Ux, Uy);
                //get the 4 neighbours cells around position (x,y)
                int x0 = Mathf.Clamp(Mathf.FloorToInt(prevParticlePos.x), 1, N - 2);
                int x1 = Mathf.Clamp(x0 + 1, 0, N - 1);
                int y0 = Mathf.Clamp(Mathf.FloorToInt(prevParticlePos.y), 1, N - 2);
                int y1 = Mathf.Clamp(y0 + 1, 0, N - 1);
                //get the deltas between x/x0 and y/y0 to weight the linear combination for the density
                float wx1 = prevParticlePos.x - x0;
                float wx0 = 1 - wx1;
                float wy1 = prevParticlePos.y - y0;
                float wy0 = 1 - wy1;
                //set back the interpolated value based on these 4 cells back at the initial grid cell (i,j)                
                s[i, j] = wx0 * (wy0 * s0[x0, y0] + wy1 * s0[x0, y1]) + wx1 * (wy0 * s0[x1, y0] + wy1 * s0[x1, y1]);
            }
        }
    }

    //trace particle back to previous position
    Vector2 TraceParticle(int i, int j, float[,] Ux, float[,] Uy)
    {
        //tracing the particle backwards in time 
        float a = (i - (dt * N) * Ux[i, j]); //new vertical position from center of cell
        float b = (j - (dt * N) * Uy[i, j]); //new horizontal position from center of cell
        // make sure the position of the particle is inside the grid; if not set it at the limit
        a = Mathf.Clamp(a,0.99f,N-0.01f);
        b = Mathf.Clamp(b,0.99f,N-0.01f);
        // return traced position
        return new Vector2(a, b);
    }

    void setBoundaryConditions(BCType bctype, float[,] x)
    {
        int iBCfactor = 1;
        int jBCfactor = 1;
        switch (bctype)
        {
            case BCType.type1:
                iBCfactor = 1;
                jBCfactor = -1;
                break;
            case BCType.type2:
                iBCfactor = -1;
                jBCfactor = 1;
                break;
            case BCType.type0:
                break;
        }
        for (int i = 1; i < N - 1; i++)
        {
            //borders along the grid
            x[i, 0] = iBCfactor * x[i, 1]; // LEFT ROW 
            x[N - 1, i] = jBCfactor * x[N - 2, i]; // TOP ROW
            x[i, N - 1] = iBCfactor * x[i, N-2]; // RIGHT ROW           
            x[0, i] = jBCfactor * x[1, i]; // BOTTOM ROW
        }
        //special case for corners
        x[0, 0] = (x[1, 0] + x[0, 1]) / 2f; //BL corner
        x[N - 1, 0] = (x[N - 2, 0] + x[N - 1, 1]) / 2f; //TL corner
        x[N - 1, N - 1] = (x[N - 2, N - 1] + x[N - 1, N - 2]) / 2f; //TR corner
        x[0, N - 1] = (x[1, N - 1] + x[0, N - 2]) / 2f; //BR corner
    }

    void dissipate(float[,] s)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                s[i, j] *= (1 / (1 + dt * dissipationRate));
            }
        }
    }
    #endregion

    //NUMERICAL METHODS

    #region numericalMethods

    float centralDifferences(float[,] x, int i, int j, int direction, float dx)
    {
        if(direction==0)
        return (x[i + 1, j] - x[i - 1, j]) / (2f * dx);
        else 
        return (x[i, j+1] - x[i, j-1]) / (2f * dx);
    }

    void GaussSeidelRelaxLinearSolver(float[,] x, float[,] x0, float leftHand, float rightHand, int limitStep = 16)
    {
        for (int step = 0; step < limitStep; step++)
        {
            for (int i = 1; i < N - 1; i++)
            {
                for (int j = 1; j < N - 1; j++)
                {
                    //linear combination of "cross" neighbours
                    x[i, j] = (x0[i, j] + leftHand * (x[i + 1, j]
                                                   + x[i, j + 1]
                                                   + x[i - 1, j]
                                                   + x[i, j - 1])) / (float)rightHand;
                }
            }
        }
    }

    float discreteVelocityDivergence(float [,] Ux, float [,] Uy, int i, int j)
    {
        return (centralDifferences(Ux, i, j, 0, N) + centralDifferences(Uy, i, j, 1, N));
    }

    #endregion

    // USER CONTROL METHODS

    #region userControls

    public void addDensitySource()
    {
        Vector2Int mousePos = renderer.getMousePosOnGrid();
        addDensity(mousePos.x, mousePos.y, 100);
    }

    public void addForce()
    {
        Vector2Int mousePos = renderer.getMousePosOnGrid();
        if (mousePos.x > 0 && mousePos.x < N && mousePos.y > 0 && mousePos.y < N)
        {
            addDensity(mousePos.x, mousePos.y, 10);
            addVelocity(mousePos.x, mousePos.y, mouseSpeed.x * mouseSpeedFactor, mouseSpeed.y * mouseSpeedFactor);
        }
    }
    void addDensity(int x, int y, float amount)
    {
        dens[x, y] += amount;
    }

    void addVelocity(int x, int y, float amountX, float amountY)
    {
        this.velX[x, y] += amountX;
        this.velY[x, y] += amountY;
    }

    void debugVelocities()
    {
        for (int i = 1; i < N - 1; i++)
        {
            for (int j = 1; j < N - 1; j++)
            {
                Vector3 initialPoint = renderer.getWorldPosFromGrid(i, j);
                Vector3 vel = initialPoint + new Vector3(1, 1, 0);
                GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);
                cube.transform.position = initialPoint + new Vector3(transform.localScale.x, transform.localScale.y, 0) * 0.5f;
                cube.transform.rotation = Quaternion.LookRotation(new Vector3(velX[i, j], velY[i, j], 0));
                cube.transform.localScale = new Vector3(0.3f, 1.5f, 0.3f);
                cube.GetComponent<Renderer>().material.color = Color.red;
                cube.transform.parent = debugVelocity.transform;
            }
        }
    }

    void clearVelocities()
    {
        for (int i = 0; i < debugVelocity.transform.childCount; i++)
        {
            GameObject.Destroy(debugVelocity.transform.GetChild(i).gameObject);
        }
    }
    #endregion
}
