using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class StableFluid : MonoBehaviour
{
    enum BCType
    {
        u_vel,
        v_vel,
        dens
    }

    //SIMULATION PARAMETERS
    [SerializeField]
    bool StartSimulation = false;
    //renderer
    [SerializeField]
    SimpleRenderer renderer;

    GameObject debugVelocity;

    //parameters

    //SOLVER PARAMETERS
    [SerializeField]
    static int N = 512;    
    [SerializeField]
    float dt = 0.2f;
    float timer = 0;

    //FLUID PARAMETERS
    [SerializeField]
    float diffRate = 0f;
    [SerializeField]
    float visc = 0.000001f;
    //2D arrays for veloctiy and density
    float[,] dens = new float[N+2, N+2];
    float[,]  dens_prv = new float[N+2, N+2];
    float[,]  u_vel = new float[N+2, N+2];
    float[,]  u_vel_prv = new float[N+2, N+2];
    float[,]  v_vel = new float[N+2, N+2];
    float[,]  v_vel_prv = new float[N+2, N+2];

    [SerializeField]
    bool debugNotDone = true;
    //handling mouse input
    Vector2 lastMousePoint;
    Vector2 mouseSpeed;
    [SerializeField]
    float mouseSpeedFactor = 30;

    //properties
    public float timeStep { get => dt; set => dt = value; }

    void Start()
    {
        //initializing the arrays
        debugVelocity = new GameObject();
        debugVelocity.transform.position = new Vector3(0, 0, 0);
        debugVelocity.name = "debugVelocities";

        addDensity(64, 64, 100);        
        renderer.renderScalarField(ref dens);
    }

    void addSourceArea(float[,] source,  float[,] x)
    {
        for (int i = 1; i < N + 1; i++)
        {
            for (int j = 1; j < N + 1; j++)
            {
                x[i, j] += dt * source[i, j];
            }
        }
    }

    void addSourcePoint(float[,] x, int a, int b, float amount)
    {
        x[a, b] += amount;
    }

    //advect the velocity/density value on the grid
    void advect( float[,] x,  float[,] x0,  float[,] u,  float[,]v)
    {
        for (int i = 1; i <N+1; i++)
        {
            for (int j = 1; j <N+1; j++)
            {
                //tracing the particle backwards in time 
                float a = ((float)i - (dt*N)*u[i, j]); //new vertical position from center of cell
                float b = ((float)j - (dt*N)*v[i, j]); //new horizontal position
                // checking if the position of the particle is inside the grid; if not set it at the limit
                a = Mathf.Max(0.5f, a);
                a = Mathf.Min(a, N + 0.5f);
                b = Mathf.Max(0.5f, b);
                b = Mathf.Min(b, N + 0.5f);
                //get the 4 neighbours cells around position (a,b)
                int a0 = Mathf.FloorToInt(a); int a1 = a0 + 1;
                int b0 = Mathf.FloorToInt(b); int b1 = b0 + 1;
                //get the deltas between a/a0 and b/b0 to weight the linear combination for the density
                float wa1 = a - a0; 
                float wa0 = 1 - wa1;
                float wb1 = b - b0; 
                float wb0 = 1 - wb1;
                //set back the interpolated value based on these 4 cells back at the initial grid cell (i,j)
                x[i, j] = wa0 * (wb0 * x0[a0, b0] + wb1 * x0[a0, b1]) + wa1 * (wb0 * x0[a1, b0] + wb1 * x0[a1, b1]);
            }
        }
    }

    //solve Poisson equation (maintain incompressibility)
    void project( float[,]u,  float [,] v,  float[,] p, float[,] div)
    {
        //compute divergence in the grid
        for (int j = 1; j < N+1; j++)
        {
            for (int i = 1; i < N+1; i++)
            {
                div[i, j] = -(0.5f) * (u[i + 1, j] - u[i - 1, j] + v[i, j + 1] - v[i, j - 1])/N;
                p[i, j] = 0;
            }            
        }
        //reset boundary conditions since we made changes in the grid
        setBC(BCType.dens,  div);
        setBC(BCType.dens,  p);
        //solve Poisson equation using GaussianSeidel solver
        GaussSeidelLinSolve( p, div, 1, 4);        
        //Hodge decomposition : any field = mass conserving + gradient
        //we want ONLY the mass conserving part (= initial field - gradient)
        for (int j = 1; j < N+1; j++)
        {
            for (int i = 1; i < N+1; i++)
            {
                u[i, j] -= 0.5f * (p[i + 1, j] - p[i - 1, j]) * N;
                v[i, j] -= 0.5f * (p[i, j+1] - p[i, j-1]) * N;
            }
        }

    }

    //diffuse density values in the grid
    void diffuse( float[,] x0,  float[,] x,float diffcoef)
    {
        float a = diffcoef * dt * N * N;
        GaussSeidelLinSolve( x,  x0, a, 1 + 4 * a);
    }

    //set boundary conditions (walls)
    void setBC(BCType bc,  float[,] x)
    {
        int iBCfactor = 1;
        int jBCfactor = 1;
        switch (bc)
        {
            case BCType.u_vel:
                iBCfactor = -1;
                jBCfactor = 1;
                break;
            case BCType.v_vel:
                iBCfactor = 1;
                jBCfactor = -1;
                break;
            case BCType.dens:
                break;
        }
        for (int i = 1; i <= N-1; i++)
        {
            //borders along the grid
            x[0, i] = iBCfactor*x[1, i]; // TOP ROW
            x[N + 1, i] = iBCfactor*x[N, i]; // BOTTOM ROW
            x[i, 0] = jBCfactor *(- x[i, 1]); // LEFT ROW 
            x[i, N + 1] = jBCfactor *(- x[i, 1]); // RIGHT ROW           
        }
        //special case for corners
        x[0, 0] = (x[1, 0] + x[0, 1]) / 2f;
        x[0, N + 1] = (x[1, N + 1] + x[0, N]) / 2f;
        x[N + 1, 0] = (x[N, 0] + x[N + 1, 1]) / 2f;
        x[N+1, N+1] = (x[N, N + 1] + x[N + 1, N]) / 2f;
    }

    //linear solver for diffuse and project steps
    void GaussSeidelLinSolve( float[,] x,  float[,] y, float numCoef, float denomCoef, int limitStep = 16)
    {
        for (int step = 0; step < limitStep; step++)
        {
            for (int i = 1; i < N+1; i++)
            {
                for (int j = 1; j < N+1; j++)
                {
                    //linear combination of "cross" neighbours
                    x[i, j] = (y[i, j] + numCoef * (x[i + 1, j]
                                                   + x[i, j + 1]
                                                   + x[i - 1, j]
                                                   + x[i, j - 1])) / (float)denomCoef;
                }
            }
        }
    }

    void densStep( float[,] dens,  float[,] dens0,  float[,] u ,  float[,] v)
    {
        addSourceArea(dens0,  dens);
        swap( dens0,  dens);
        diffuse( dens0,  dens, diffRate);
        setBC(BCType.dens, dens);
        swap( dens0, dens);        
        advect( dens, dens0, u, v);
        setBC(BCType.dens,  dens);

    }

    void velStep( float[,] u, float[,] u0,   float[,] v,  float[,] v0)
    {
        //add force
        addSourceArea(u0,  u); addSourceArea(v0,  v);
        //diffusion + reset bc
        swap( u0,  u); 
        diffuse(u0, u, visc); setBC(BCType.u_vel, u);
        swap( v0,  v); 
        diffuse(v0, v, visc); setBC(BCType.v_vel, v);
        //projection
        project(u0, v0, u, v);
        setBC(BCType.u_vel, u);
        setBC(BCType.v_vel, v);
        swap( u0,  u); swap( v0,  v);
        //advection + reset bc
        advect( u, u0, u0, v0); setBC(BCType.u_vel, u);
        advect( v, v0, u0, v0); setBC(BCType.v_vel, v);
        //projection
        project(u, v, u0, v0);
        setBC(BCType.u_vel, u);
        setBC(BCType.v_vel, v);
    }

    void fluidStep()
    {
        velStep(u_vel,u_vel_prv,v_vel,v_vel_prv);
        densStep(dens,dens_prv,u_vel,v_vel);
    }

    void Update()
    {
        if (StartSimulation)
        {
            timer += Time.deltaTime;
            if (timer > timeStep)
            {
                fluidStep();
                timer -= timeStep;
                renderer.renderScalarField(ref dens);
            }
        }

        if (Input.GetKeyDown(KeyCode.V))
        {
            debugVelocities();            
        }

        if (Input.GetKeyDown(KeyCode.U))
        {
            debugOldVelocities();
        }

        if (Input.GetKeyDown(KeyCode.Escape))
        {
            clearVelocities();
        }
    }

    void swap( float[,] v1,  float[,] v2)
    {
        float[,] temp = v2;
        v2 = v1;
        v1 = temp;
    }


    public void addDensity()
    {
        Vector2Int mousePos = renderer.getMousePosOnGrid();
        addSourcePoint(dens_prv, mousePos.x, mousePos.y, 100);
    }

    void addVelocity(int x, int y, float amountX, float amountY)
    {        
        u_vel_prv[x,y] += amountX *dt;
        v_vel_prv[x,y] += amountY*dt;
    }

    void addDensity(int x, int y, float amount)
    {        
        dens_prv[x,y] += amount*dt;
    }


    public void addVelocity()
    {        
        Vector2Int mousePos = renderer.getMousePosOnGrid();
        addSourcePoint(u_vel_prv, mousePos.x, mousePos.y, mouseSpeed.x * mouseSpeedFactor );
        addSourcePoint(v_vel_prv, mousePos.x, mousePos.y, mouseSpeed.y * mouseSpeedFactor );
        addSourcePoint(dens_prv, mousePos.x, mousePos.y, 20);
    }

    void debugVelocities()
    {
        for (int i = 1; i < N + 1; i++)
        {
            for (int j = 1; j < N + 1; j++)
            {
                Vector3 initialPoint = renderer.getWorldPosFromGrid(i, j);
                Vector3 vel = initialPoint + new Vector3(1, 1, 0);
                GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);
                cube.transform.position = initialPoint + new Vector3(transform.localScale.x, transform.localScale.y, 0) * 0.5f;
                cube.transform.rotation = Quaternion.LookRotation(new Vector3(u_vel[i, j], v_vel[i, j], 0));
                cube.transform.localScale = new Vector3(0.3f, 1.5f, 0.3f);
                cube.GetComponent<Renderer>().material.color = Color.red;
                cube.transform.parent = debugVelocity.transform;
            }
        }
    }

    void debugOldVelocities()
    {
        for (int i = 1; i < N + 1; i++)
        {
            for (int j = 1; j < N + 1; j++)
            {
                Vector3 initialPoint = renderer.getWorldPosFromGrid(i, j);
                Vector3 vel = initialPoint + new Vector3(1, 1, 0);
                GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);
                cube.transform.position = initialPoint + new Vector3(transform.localScale.x, transform.localScale.y, 0) * 0.5f;
                cube.transform.rotation = Quaternion.LookRotation(new Vector3(u_vel_prv[i, j], v_vel_prv[i, j], 0));
                cube.transform.localScale = new Vector3(0.3f, 1.5f, 0.3f);
                cube.GetComponent<Renderer>().material.color = Color.blue;
                cube.transform.parent = debugVelocity.transform;
            }
        }
    }

    void printSum(float[,] x)
    {
        float sum = 0;
        for(int i=1;i<N+1;i++)
        {
            for (int j = 1; j < N + 1; j++)
            {
                sum += x[i, j];
            }
        }
        print("sum = " + sum);
    }

    void clearVelocities()
    {
        for (int i = 0; i < debugVelocity.transform.childCount; i++)
        {
            GameObject.Destroy(debugVelocity.transform.GetChild(i).gameObject);
        }
    }
}
