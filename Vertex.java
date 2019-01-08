
public class Vertex {
	protected StringBuffer _strUniqueID, //a unique id identifying vertex
	 _strData; //data associated with vertex
	 protected int _nX,_nY; 
	 
	 public Vertex(StringBuffer _strUniqueID,StringBuffer _strData,int _nX,int _nY)
	 {
		 this._nX=_nX;
		 this._nY=_nY;
		 this._strUniqueID=_strUniqueID;
		 this._strData=_strData;
	 }
	 public StringBuffer getUniqueID( )
	 {
		 return _strUniqueID;
		  }

		  public StringBuffer getData( )
		  {
		  return _strData;
		  }
		  public int getX( )
		  {
		  return _nX;
		  }
		  public int getY( )
		  {
		  return _nY;
		  }
}
