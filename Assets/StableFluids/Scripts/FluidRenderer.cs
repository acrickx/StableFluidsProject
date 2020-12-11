using UnityEngine;
using UnityEngine.UI;

public class FluidRenderer : MonoBehaviour
{
    [SerializeField]
    private int resolution = 128;
    [SerializeField]
    private RawImage renderImage;
    //[SerializeField]
    //private Image[] corners;
    [SerializeField]
    private Image BL;
    [SerializeField]
    private Camera cam;
    private Texture2D renderTexture;
    [SerializeField]
    [Range(1,32)]
    private int currentZ=1;

    void Awake()
    {                 
        renderTexture = new Texture2D(resolution, resolution);
        //renderTexture.filterMode = FilterMode.Point;
        renderImage.texture = renderTexture;        
    }

    public void renderScalarField(ref float[,] scalarField)
    {
        for (int i = 0; i < resolution; i++)
        {
            for (int j = 0; j < resolution; j++)
            {
                renderTexture.SetPixel(i, j, new Color(scalarField[i, j], scalarField[i, j], scalarField[i, j]));
            }
        }
        renderTexture.Apply();
    }

    public void renderScalarField3D(ref float[,,] scalarField)
    {
        for (int i = 0; i < resolution; i++)
        {
            for (int j = 0; j < resolution; j++)
            {
                renderTexture.SetPixel(i, j, new Color(scalarField[i, j,currentZ], scalarField[i, j,currentZ], scalarField[i, j,currentZ]));
            }
        }
        renderTexture.Apply();
    }

    public Vector2Int getMousePosOnGrid()
    {
        float minX, minY;
        minX = cam.WorldToScreenPoint(BL.transform.position).x;
        minY = cam.WorldToScreenPoint(BL.transform.position).y;
        Vector2Int pos = new Vector2Int(Mathf.RoundToInt((Input.mousePosition.x-minX)/renderImage.transform.localScale.x), Mathf.RoundToInt((Input.mousePosition.y - minY) / renderImage.transform.localScale.y));
        return pos;    
    }

    public Vector3Int getMousePosOnGrid3D()
    {
        float minX, minY;
        minX = cam.WorldToScreenPoint(BL.transform.position).x;
        minY = cam.WorldToScreenPoint(BL.transform.position).y;
        Vector3Int pos = new Vector3Int(Mathf.RoundToInt((Input.mousePosition.x - minX) / renderImage.transform.localScale.x), Mathf.RoundToInt((Input.mousePosition.y - minY) / renderImage.transform.localScale.y), currentZ);
        return pos;
    }

    public Vector3 getWorldPosFromGrid(float i, float j)
    {
        return cam.ScreenToWorldPoint(cam.WorldToScreenPoint(BL.transform.position) + new Vector3(i * renderImage.transform.localScale.x, j * renderImage.transform.localScale.y, 0));       
    }

}
