using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LiquidTextureAnimator : MonoBehaviour
{
    enum BCType
    {
        vertical,
        horizontal,
        pressure
    }

    //SIMULATION PARAMETERS
    [SerializeField]
    bool StartSimulation = false;
    //renderer
    [SerializeField]
    RawImage image;
    [SerializeField]
    Texture2D originalTexture;

    //SOLVER PARAMETERS
    [SerializeField]
    int N = 512;
    float dx = 0;
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

    //RENDERING
    [SerializeField]
    private Image BL;
    [SerializeField]
    private Camera cam;

    //2D arrays for veloctiy and dens
    float[,] xCoord;
    float[,] xCoordPrev;
    float[,] yCoord;
    float[,] yCoordPrev;
    float[,] velY;
    float[,] velYPrev;
    float[,] velX;
    float[,] velXPrev;

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
        xCoord = new float[N, N];
        xCoordPrev = new float[N, N];
        yCoord = new float[N, N];
        yCoordPrev = new float[N, N];
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                xCoord[i, j] = i / (float)N;
                yCoord[i, j] = j / (float)N;
            }
        }
        renderLiquidTexture();

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
        diffusionStep(velXPrev, velX, viscosity);
        setBoundaryConditions(BCType.horizontal, velX);
        diffusionStep(velYPrev, velY, viscosity);
        setBoundaryConditions(BCType.vertical, velY);
        //projection step
        projectionStep(velXPrev, velYPrev, velY, velX);
        //Advection step
        advectionStepFE(velX, velXPrev, velYPrev, velXPrev);
        setBoundaryConditions(BCType.horizontal, velX);
        advectionStepFE(velY, velYPrev, velYPrev, velXPrev);
        setBoundaryConditions(BCType.vertical, velY);
        //projection step -> solving Poisson equation
        projectionStep(velY, velX, velYPrev, velXPrev);
    }

    void Sstep()
    {
        //diffusion step
        diffusionStep(xCoordPrev, xCoord, diffusionRate);
        diffusionStep(yCoordPrev, yCoord, diffusionRate);
        //advection step
        advectionStepFE(xCoord, xCoordPrev, velY, velX);
        advectionStepFE(yCoord, yCoordPrev, velY, velX);
        setBoundaryConditions(BCType.pressure, xCoord);
        setBoundaryConditions(BCType.pressure, yCoord);
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
                renderLiquidTexture();
            }
        }
        mouseSpeed = new Vector2(Input.mousePosition.x, Input.mousePosition.y) - lastMousePoint;
        lastMousePoint = Input.mousePosition;
        print("FPS : " + 1 / Time.deltaTime);
    }

    #endregion

    //SOLVER METHODS

    #region solverMethods

    void diffusionStep(float[,] x, float[,] x0, float diffusion)
    {
        float diffusionRate = dt * diffusion / (dx * dx);
        NumericalMethods.GaussSeidelRelaxLinearSolver(x, x0, diffusionRate, 1 + 4 * diffusionRate, N);        
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
                div[i, j] = NumericalMethods.divergence(Ux, Uy, i, j, N) * (-dx * dx);
                //setting initial to pressure for linear solve
                pressure[i, j] = 0;
            }
        }
        //reset bondary conditions since we made changes in the grid
        setBoundaryConditions(BCType.pressure, div); setBoundaryConditions(BCType.pressure, pressure);
        //solve Poisson equation using linear solver to find correct pressure values
        NumericalMethods.GaussSeidelRelaxLinearSolver(pressure, div, 1, 4, N);        
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
                Vector2 prevParticlePos = NumericalMethods.TraceParticleFE(i, j, Ux, Uy, dt, N);
                //interpolate from adjacent nodes the scalar field value and set it in s
                s[i, j] = NumericalMethods.bilinearInterpolation(s0, prevParticlePos, N);
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
            x[i, N - 1] = iBCfactor * x[i, N - 2]; // RIGHT ROW           
            x[0, i] = jBCfactor * x[1, i]; // BOTTOM ROW
        }
        //special case for corners
        x[0, 0] = (x[1, 0] + x[0, 1]) / 2f; //BL corner
        x[N - 1, 0] = (x[N - 2, 0] + x[N - 1, 1]) / 2f; //TL corner
        x[N - 1, N - 1] = (x[N - 2, N - 1] + x[N - 1, N - 2]) / 2f; //TR corner
        x[0, N - 1] = (x[1, N - 1] + x[0, N - 2]) / 2f; //BR corner
    }

    #endregion
    // USER CONTROL METHODS

    #region userControls

    public void addForce()
    {
        Vector2Int mousePos = getMousePosOnGrid();
        if (mousePos.x > 0 && mousePos.x < N && mousePos.y > 0 && mousePos.y < N)
        {            
            addVelocity(mousePos.x, mousePos.y, mouseSpeed.y * mouseSpeedFactor, mouseSpeed.x * mouseSpeedFactor);
        }
    }

    void addVelocity(int x, int y, float amountX, float amountY)
    {
        this.velX[x, y] += amountX;
        this.velY[x, y] += amountY;
    }

    public Vector3 getWorldPosFromGrid(float i, float j)
    {
        return cam.ScreenToWorldPoint(cam.WorldToScreenPoint(BL.transform.position) + new Vector3(i * image.transform.localScale.x, j * image.transform.localScale.y, 0));
    }

    void renderLiquidTexture()
    {
        Texture2D updatedTexture = new Texture2D(N, N);
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                updatedTexture.SetPixel(i, j, originalTexture.GetPixel(Mathf.RoundToInt(xCoord[i, j] * N), Mathf.RoundToInt(yCoord[i, j] * N)));
            }
        }
        updatedTexture.Apply();
        image.texture = updatedTexture;
    }

    public Vector2Int getMousePosOnGrid()
    {
        float minX, minY;
        minX = cam.WorldToScreenPoint(BL.transform.position).x;
        minY = cam.WorldToScreenPoint(BL.transform.position).y;
        Vector2Int pos = new Vector2Int(Mathf.RoundToInt((Input.mousePosition.x - minX) / image.transform.localScale.x), Mathf.RoundToInt((Input.mousePosition.y - minY) / image.transform.localScale.y));
        return pos;
    }

    #endregion
}
