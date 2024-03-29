﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class StableFluidSolver : MonoBehaviour
{
    enum BCType
    {
        vertical,
        horizontal,
        pressure
    }

    enum advectionTechnique
    {
        FE,
        RK2,
        RK3
    }

    enum linearSolver
    {
        Jacobi,
        GaussSeidel
    }

    //SIMULATION PARAMETERS
    [SerializeField]
    bool StartSimulation = false;
    //renderer
    [SerializeField]
    FluidRenderer renderer;

    //SOLVER PARAMETERS
    [SerializeField]
    int N = 128;
    float dx =0;
    [SerializeField]
    float dt = 0.2f;
    float timer = 0;
    [SerializeField]
    advectionTechnique technique;
    [SerializeField]
    linearSolver solver;

    //FLUID PARAMETERS
    [SerializeField]
    float dissipation = 0.2f;
    [SerializeField]
    float diffusion = 0f;
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

        dx = 1 / (float)N;

        debugVelocity = new GameObject();
        debugVelocity.name = "debugVelocity";
    }

    //MAIN LOOP
    #region mainLoop

    void Vstep()
    {
        //difusion step
        diffusionStep(velXPrev, velX, viscosity); 
        setBoundaryConditions(BCType.horizontal, velX);

        diffusionStep(velYPrev, velY, viscosity); 
        setBoundaryConditions(BCType.vertical, velY);
        
        //projection step -> keep divergence-free velocity
        projectionStep(velXPrev, velYPrev, velX, velY);
        setBoundaryConditions(BCType.vertical, velY);
        setBoundaryConditions(BCType.horizontal, velX);

        //self-advection step -> different numerical techniques        
        if (technique == advectionTechnique.RK3)
        {
            advectionStepRK3(velX, velXPrev, velXPrev, velYPrev);
            setBoundaryConditions(BCType.horizontal, velX);            

            advectionStepRK3(velY, velYPrev, velXPrev, velYPrev);
            setBoundaryConditions(BCType.vertical, velY);            
        }
        else if (technique == advectionTechnique.RK2)
        {
            advectionStepRK2(velX, velXPrev, velXPrev, velYPrev);
            setBoundaryConditions(BCType.horizontal, velX);

            advectionStepRK2(velY, velYPrev, velXPrev, velYPrev);
            setBoundaryConditions(BCType.vertical, velY);            
        }
        else
        {
            advectionStepFE(velX, velXPrev, velXPrev, velYPrev);
            setBoundaryConditions(BCType.horizontal, velX);

            advectionStepFE(velY, velYPrev, velXPrev, velYPrev);
            setBoundaryConditions(BCType.vertical, velY);            
        }

        //projection step -> keep divergence-free velocity
        projectionStep(velX, velY, velXPrev, velYPrev);
        setBoundaryConditions(BCType.horizontal, velY);
        setBoundaryConditions(BCType.vertical, velX);
    }

    void Sstep()
    {
        //diffusion step
        diffusionStep(densPrev, dens, diffusion);

        //advection step
        if (technique == advectionTechnique.FE)
        {
            advectionStepFE(dens, densPrev, velX, velY);
        }
        else if (technique == advectionTechnique.RK2)
        {
            advectionStepRK2(dens, densPrev, velX, velY);
        }
        else
        {
            advectionStepRK3(dens, densPrev, velX, velY);         
        }

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
        print("FPS : " + 1 / Time.deltaTime);
    }

    #endregion

    //SOLVER METHODS

    #region solverMethods

    void diffusionStep(float[,] x, float[,] x0, float diffusion)
    {
        float diffusionRate = dt * diffusion / (dx*dx);
        if (solver == linearSolver.Jacobi)
        {
            NumericalMethods.JacobiLinearSolver(x, x0, diffusionRate, 1 + 4 * diffusionRate, N);
        }
        else
        {
            NumericalMethods.GaussSeidelRelaxLinearSolver(x, x0, diffusionRate, 1 + 4 * diffusionRate, N);
        }
    }

    void projectionStep(float[,] Ux, float[,] Uy, float[,] pressure, float[,] div)
    {
        //compute divergence in the grid
        for (int i = 1; i < N - 1; i++)
        {
            for (int j = 1; j < N - 1; j++)
            {
                //computing divergence which will be the right hand of the equation for pressure linear system
                //(-1*dx*dx) because we isolate pressure[i,j] in the equation for the solver
                div[i, j] = NumericalMethods.divergence(Ux, Uy, i, j,N) * (-dx*dx);
                //setting initial to pressure for linear solve
                pressure[i, j] = 0;
            }
        }
        //reset bondary conditions since we made changes in the grid
        setBoundaryConditions(BCType.pressure, div); setBoundaryConditions(BCType.pressure, pressure);
        //solve Poisson equation using linear solver to find correct pressure values
        if (solver == linearSolver.Jacobi)
        {
            NumericalMethods.JacobiLinearSolver(pressure, div, 1, 4, N);
        }
        else
        {
            NumericalMethods.GaussSeidelRelaxLinearSolver(pressure, div, 1, 4, N);
        }
        //subtract pressure gradient
        project(Ux, Uy, pressure);
    }

    //Hodge decomposition : any field = mass conserving + gradient
    //we want ONLY the mass conserving part (= initial field - gradient)
    void project(float[,] Ux, float[,] Uy, float[,] pressure)
    {
        for (int j = 1; j < N - 1; j++)
        {
            for (int i = 1; i < N - 1; i++)
            {
                Vector2 pressureGrad = NumericalMethods.grad(pressure, i, j, dx);
                Ux[i, j] -= pressureGrad.x;
                Uy[i, j] -= pressureGrad.y;
            }
        }
    }

    void advectionStepFE(float[,] s, float[,] s0, float[,] Ux, float[,] Uy)
    {
        for (int i = 1; i < N - 1; i++)
        {
            for (int j = 1; j < N - 1; j++)
            {
                //trace "imaginary" particle from the grid position back in time to previous position
                Vector2 prevParticlePos = NumericalMethods.TraceParticleFE(i, j, Ux, Uy,dt,N);                
                //interpolate from adjacent nodes the scalar field value and set it in s
                s[i, j] = NumericalMethods.bilinearInterpolation(s0, prevParticlePos,N);                
            }
        }
    }

    void advectionStepRK2(float[,] s, float[,] s0, float[,] Ux, float[,] Uy)
    {
        for (int i = 1; i < N - 1; i++)
        {
            for (int j = 1; j < N - 1; j++)
            {
                //trace "imaginary" particle from the grid position back in time to previous position
                Vector2 prevParticlePos = NumericalMethods.TraceParticleRK2(i, j, Ux, Uy, dt, N);
                //interpolate from adjacent nodes the scalar field value and set it in s
                s[i, j] = NumericalMethods.bilinearInterpolation(s0, prevParticlePos, N);
            }
        }
    }

    void advectionStepRK3(float[,] s, float[,] s0, float[,] Ux, float[,] Uy)
    {
        for (int i = 1; i < N - 1; i++)
        {
            for (int j = 1; j < N - 1; j++)
            {
                //trace "imaginary" particle from the grid position back in time to previous position
                Vector2 prevParticlePos = NumericalMethods.TraceParticleRK3(i, j, Ux, Uy, dt, N);
                //interpolate from adjacent nodes the scalar field value and set it in s
                s[i, j] = NumericalMethods.bilinearInterpolation(s0, prevParticlePos, N);
            }
        }
    }

    void setBoundaryConditions(BCType bctype, float[,] x)
    {
        int iBCfactor = 1;
        int jBCfactor = 1;
        switch (bctype)
        {
            case BCType.vertical:
                iBCfactor = -1;
                jBCfactor = 1;
                break;
            case BCType.horizontal:
                iBCfactor = 1;
                jBCfactor = -1;
                break;
            case BCType.pressure:
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
                s[i, j] *= (1 / (1 + dt * dissipation));
            }
        }
    }
    #endregion

    // USER CONTROL METHODS

    #region userControls

    public void addDensitySource()
    {
        Vector2Int mousePos = renderer.getMousePosOnGrid();
        addDensity(mousePos.x, mousePos.y, 50);
    }

    public void addForce()
    {
        Vector2Int mousePos = renderer.getMousePosOnGrid();
        if (mousePos.x > 0 && mousePos.x < N && mousePos.y > 0 && mousePos.y < N)
        {
            addDensity(mousePos.x, mousePos.y, 10);
            addVelocity(mousePos.x, mousePos.y, mouseSpeed.y * mouseSpeedFactor, mouseSpeed.x * mouseSpeedFactor);
        }
    }
    void addDensity(int x, int y, float amount)
    {
        dens[x, y] += amount;
    }

    void addVelocity(int x, int y, float amountX, float amountY)
    {
        this.velY[x, y] += amountX;
        this.velX[x, y] += amountY;
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
                cube.transform.rotation = Quaternion.LookRotation(new Vector3(velY[i, j], velX[i, j], 0));
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
