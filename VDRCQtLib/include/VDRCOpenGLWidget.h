#pragma once

#include <QtOpenGL>
#include <gl/GLU.h>

#include "constForVDRCOpenGLWidget.h"

#include "rg_Point3D.h"
#include "rg_TMatrix3D.h"
#include "BallGeneratorCore.h"
#include "VEdgeCore.h"
#include "VVertexCore.h"
#include "VFaceCore.h"
#include "Ellipsoid3D.h"

#include "Color3f.h"
#include <list>
#include <vector>
#include <array>
#include <map>
#include <set>

using namespace std;
using namespace V;
using namespace CoreTier;

class VDRCOpenGLWidget : public QOpenGLWidget
{
	Q_OBJECT

public:
	//rg_TMatrix3D m_transform;

	rg_Point3D localOrigin;
	rg_Point3D localX;
	rg_Point3D localY;
	rg_Point3D localZ;

	GLfloat positionLight[2][4];
	GLfloat specularLight[4];
	GLfloat diffuseLight[4];
	GLfloat ambientLight[4];

	GLfloat ambientMaterial[4];
	GLfloat specularMaterial[4];
	GLfloat emittedMaterial[4];
	GLfloat shininess;

	GLUquadricObj* qObj;

	double coefficientsForOthogornalPlaneFromCamera[4];

	float zoomFactor;

	QPoint lastPos;

	//For picking
	PICK_MODE m_pickMode;
	GLuint m_buff[SELECTION_BUFFER_SIZE];
	GLint   m_hits;
	map<int, VVertexCore*> m_mapForVVertexID;

	VVertexCore* m_selectedVtx;
	bool m_bDisableRotation;

	int m_eyeDistance = 1000;

public:
	explicit VDRCOpenGLWidget(QWidget *parent = 0);
	virtual ~VDRCOpenGLWidget();

	void initializeGL();
	void paintGL();
	void resizeGL(int width, int height);
	void perspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar);

	virtual void draw() = 0;

	QSize minimumSizeHint() const;
	QSize sizeHint() const;
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent* event);
	

	inline void set_eye_distance(const float& distance) { m_eyeDistance = distance; }
	inline void set_local_origin(const float& x, const float& y, const float& z) { localOrigin.setPoint(x, y, z); }

	bool		initialize_eye_position();
	bool    rotate_eye_position(float angleX, float angleY, float angleZ);
	//bool		rotate_eye_position_at_local_cener_point(float angleX, float angleY, float angleZ, rg_TMatrix3D& localCenterPt);
	void		set_eye_direction(rg_Point3D eyeDirection);

	bool    zoom(const float fScale);

	void    get_rotation(const rg_TMatrix3D& rotationMatrix, float*& rotationFactors) const;

	//Drawing voronoi elements
	void draw_generators(const list<BallGeneratorCore*>& generators) const;
	void draw_voronoi_vertex(const list<VVertexCore*>& VVertices, const float& ballRadius, const Color3f& color, const float& A = 1.0);
	void draw_voronoi_edges(const list<VEdgeCore*>& VEdges, const float& thickness, const Color3f& color, const float& A = 1.0, const bool& isStipple = false) const;
		void draw_voronoi_edge(const VEdgeCore* VEdge, const float& thickness, const Color3f& color, const float& A = 1.0, const bool& isStipple = false) const;
	void draw_voronoi_faces(const list<VFaceCore*>& VFaces, const Color3f& color, const float& A = 1.0) const;

	//Drawing geometric elements
	void draw_sphere(const rg_Point3D& center, const float& radius, const Color3f& color, const float& A = 1.0, const int& elementID = -1) const;
	void draw_ellipsoid(const Ellipsoid3D& ellipsoid, const Color3f& color, const float& A = 1.0, const int& elementID = -1) const;
	void draw_point(rg_Point3D& pt, const float& ptSize, const Color3f& color, const float& A = 1.0, const int& elementID = -1) const;
	void draw_line(const rg_Point3D& pt1, const rg_Point3D& pt2, const float& width, const Color3f& color, const float& A = 1.0) const;
	void draw_line_stipple(const rg_Point3D& pt1, const rg_Point3D& pt2, const float& thickness, const Color3f& color, const float& A = 1.0) const;
	void draw_face(const list<rg_Point3D>& points, const Color3f& color, const float& A = 1.0) const;
	void draw_triangle(const array<rg_Point3D, 3>& points, const Color3f& color, const float& A = 1.0) const;
	void draw_octagonal_cone(const rg_Point3D& base, const rg_Point3D& tip, const float& radius, const Color3f& color, const float& A = 1.0) const;

	//For picking
	void gl_select(int x, int y);
	bool find_selected_elements();
	VVertexCore* find_closest_VVertex(const list<int>& elementIDList);
	void enlist_VVertices(set<VVertexCore*>& vertices);
};