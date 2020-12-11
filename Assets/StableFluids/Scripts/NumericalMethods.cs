using System.Collections;
using System.Collections.Generic;
using UnityEngine;

static class NumericalMethods
{

    //2D NUMERICAL METHODS
    #region 2DMethods

    public static float centralDifferences(float[,] x, int i, int j, int direction, float dx)
    {
        if (direction == 0)
            return (x[i + 1, j] - x[i - 1, j]) / (2f * dx);
        else
            return (x[i, j + 1] - x[i, j - 1]) / (2f * dx);
    }

    public static Vector2 grad(float[,] x, int i, int j, float dx)
    {
        return new Vector2(centralDifferences(x, i, j, 1,dx), centralDifferences(x,i,j,0,dx));        
    }

    public static float bilinearInterpolation(float[,] s, Vector2 pos, int N)
    {
        int minX = Mathf.Clamp(Mathf.FloorToInt(pos.x), 1, N - 2);
        int maxX = Mathf.Clamp(minX + 1, 0, N - 1);
        int minY = Mathf.Clamp(Mathf.FloorToInt(pos.y), 1, N - 2);
        int maxY = Mathf.Clamp(minY + 1, 0, N - 1);
        float alpha = pos.x - minX;
        float beta = pos.y - minY;
        return (1 - alpha) * ((1 - beta) * s[minX, minY] + beta * s[minX, maxY]) + alpha * ((1 - beta) * s[maxX, minY] + beta * s[maxX, maxY]);
    }

    //Linear solver from Jos Stam article Real Time Fluid Dynamics from Games
    public static void GaussSeidelRelaxLinearSolver(float[,] x, float[,] x0, float numCoef, float denumCoef, int N, int limitStep = 16)
    {
        for (int step = 0; step < limitStep; step++)
        {
            for (int i = 1; i < N - 1; i++)
            {
                for (int j = 1; j < N - 1; j++)
                {
                    x[i, j] = (x0[i, j] + numCoef * (x[i + 1, j]
                                                   + x[i, j + 1]
                                                   + x[i - 1, j]
                                                   + x[i, j - 1])) / (float)denumCoef;
                }
            }
        }
    }

    public static void JacobiLinearSolver(float[,] x, float[,] x0, float numCoef, float denumCoef, int N, int limitStep = 16)
    {
        float[,] aux = new float[N,N];
        for (int step = 0; step < limitStep; step++)
        {
            for (int i = 1; i < N - 1; i++)
            {
                for (int j = 1; j < N - 1; j++)
                {
                    aux[i,j] = (x0[i, j] + numCoef * (x[i + 1, j]
                                                   + x[i, j + 1]
                                                   + x[i - 1, j]
                                                   + x[i, j - 1])) / (float)denumCoef;
                }
            }
            for (int i = 1; i < N - 1; i++)
            {
                for (int j = 1; j < N - 1; j++)
                {
                    x[i, j] = aux[i, j];
                }
            }
        }
    }


    public static float divergence(float[,] Ux, float[,] Uy, int i, int j, int N)
    {
        return (NumericalMethods.centralDifferences(Ux, i, j, 1, 1/(float)N) + NumericalMethods.centralDifferences(Uy, i, j, 0, 1/(float)N));
    }

    //TRACE PARTICLE METHODS

    //Forward Euler
    public static Vector2 TraceParticleFE(int i, int j, float[,] Ux, float[,] Uy, float dt, int N)
    {
        //tracing the particle backwards in time 
        float j1 = (j - dt * Ux[i, j]); //new horizontal position from center of cell
        float i1 = (i - dt * Uy[i, j]); //new vertical position from center of cell
        // make sure the position of the particle is inside the grid; if not set it at the limit
        i1 = Mathf.Clamp(i1, 0.99f, N - 0.01f);
        j1 = Mathf.Clamp(j1, 0.99f, N - 0.01f);
        // return traced position
        return new Vector2(i1, j1);
    }

    //Runge-Kutta 2nd order
    public static Vector2 TraceParticleRK2(int i, int j, float[,] Ux, float[,] Uy, float dt, int N)
    {
        //tracing the particle backwards in time 
        float jmid = (j - 0.5f * dt * Ux[i, j]); //halfway horizontal position from center of cell
        float imid = (i - 0.5f * dt * Uy[i, j]); //halfway vertical position from center of cell
        Vector2 halfwayPos = new Vector2(imid, jmid);
        float j1 = j - dt * bilinearInterpolation(Ux, halfwayPos, N);
        float i1 = i - dt * bilinearInterpolation(Uy, halfwayPos, N);
        // make sure the position of the particle is inside the grid; if not set it at the limit
        i1 = Mathf.Clamp(i1, 0.99f, N - 0.01f);
        j1 = Mathf.Clamp(j1, 0.99f, N - 0.01f);
        // return traced position
        return new Vector2(i1, j1);
    }

 
    //Runge-Kutta 3rd order based on Ralston article from 1962
    public static Vector2 TraceParticleRK3(int i, int j, float[,] Ux, float[,] Uy, float dt, int N)
    {
        //tracing the particle backwards in time 
        Vector2 k1 = new Vector2(Uy[i, j], Ux[i, j]); //halfway vertical position from center of cell   
        float imid = (i - 0.5f * (dt * k1.x));
        float jmid = (j - 0.5f * (dt * k1.y));
        Vector2 halfwayPos = new Vector2(imid, jmid);
        Vector2 k2 = new Vector2(bilinearInterpolation(Uy, halfwayPos, N), bilinearInterpolation(Ux, halfwayPos, N));
        float i3quarter = (i - 0.75f * dt * k2.x);
        float j3quarter = (j - 0.75f * dt * k2.y);
        Vector2 threeQuarterPos = new Vector2(i3quarter, j3quarter);
        Vector2 k3 = new Vector2(bilinearInterpolation(Uy, threeQuarterPos, N), bilinearInterpolation(Ux, threeQuarterPos, N));
        float i1 = i - ((2f / 9f) * dt * k1.x + (3f / 9f) * dt * k2.x + (4f / 9f) * dt * k3.x);
        float j1 = j - ((2f / 9f) * dt * k1.y + (3f / 9f) * dt * k2.y + (4f / 9f) * dt * k3.y);
        // make sure the position of the particle is inside the grid; if not set it at the limit
        i1 = Mathf.Clamp(i1, 0.99f, N - 0.01f);
        j1 = Mathf.Clamp(j1, 0.99f, N - 0.01f);
        return new Vector2(i1, j1);
    }

    #endregion

    //3D NUMERICAL METHODS

    #region 3DMethods

    public static float discreteVelocityDivergence3D(float[,,] Ux, float[,,] Uy, float[,,] Uz, int i, int j, int k, int N)
    {
        return (centralDifferences3D(Ux, i, j, k, 1, N) + centralDifferences3D(Uy, i, j, k, 0, N)) + centralDifferences3D(Uz, i, j, k, 2, N);
    }

    public static float centralDifferences3D(float[,,] x, int i, int j, int k, int direction, float dx)
    {
        if (direction == 0)
            return (x[i + 1, j, k] - x[i - 1, j, k]) / (2f * dx);
        else if (direction == 1)
            return (x[i, j + 1, k] - x[i, j - 1, k]) / (2f * dx);
        else
            return (x[i, j, k + 1] - x[i, j - 1, k - 1]) / (2f * dx);
    }

    public static void GaussSeidelRelaxLinearSolver3D(float[,,] x, float[,,] x0, float leftHand, float rightHand, int N, int limitStep = 16)
    {
        for (int step = 0; step < limitStep; step++)
        {
            for (int k = 1; k < N - 1; k++)
            {
                for (int i = 1; i < N - 1; i++)
                {
                    for (int j = 1; j < N - 1; j++)
                    {
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

    public static float trilinearInterpolation(float[,,] s0, Vector3 pos, int N)
    {
        int minX = Mathf.Clamp(Mathf.FloorToInt(pos.x), 1, N - 2);
        int maxX = Mathf.Clamp(minX + 1, 0, N - 1);
        int minY = Mathf.Clamp(Mathf.FloorToInt(pos.y), 1, N - 2);
        int maxY = Mathf.Clamp(minY + 1, 0, N - 1);
        int minZ = Mathf.Clamp(Mathf.FloorToInt(pos.z), 1, N - 2);
        int maxZ = Mathf.Clamp(minZ + 1, 0, N - 1);        
        float alpha = pos.x - minX;
        float beta = pos.y - minY;        
        float gamma = pos.z - minZ;              
        return (1-alpha)* ((1-beta)* ((1-gamma)* s0[minX, minY, minZ] + gamma* (s0[minX, minY, maxZ])) + beta* ((1-gamma)* s0[minX, maxY, minZ] + gamma* (s0[minX, maxY, maxZ])))
                                   + alpha* ((1-beta)* ((1-gamma)* s0[maxX, minY, minZ] + gamma* (s0[maxX, minY, maxZ])) + beta* ((1-gamma)* s0[maxX, maxY, minZ] + gamma* (s0[maxX, maxY, maxZ])));         
    }


    //trace particle back to previous position
    public static Vector3 TraceParticle3D(int i, int j, int k, float[,,] Ux, float[,,] Uy, float[,,] Uz, float dt, int N)
    {
        //tracing the particle backwards in time 
        float a = (i - (dt * N) * Uy[i, j, k]); //new vertical position from center of cell
        float b = (j - (dt * N) * Ux[i, j, k]); //new horizontal position
        float c = (k - (dt * N) * Uz[i, j, k]); //new depth position
        // make sure the position of the particle is inside the grid; if not set it at the limit
        a = Mathf.Clamp(a, 0.99f, N - 0.01f);
        b = Mathf.Clamp(b, 0.99f, N - 0.01f);
        c = Mathf.Clamp(c, 0.99f, N - 0.01f);
        // return traced position
        return new Vector3(a, b, c);
    }
    #endregion
}
