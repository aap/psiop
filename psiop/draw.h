struct VertCol {
	Vec4 pos;
	uint8_t r, g, b, a;
};

void initDraw(uint16_t *fb, uint32_t *zb, int w, int h);
void clear(int zbuf, uint16_t col);
void updatemat(void);
void drawPoint(VertCol v);
void drawLine(VertCol v1, VertCol v2);
void drawTriangle(VertCol v1, VertCol v2, VertCol v3);
void setVertices(VertCol *scrverts, int numVerts);
void drawIndexed(int prim, uint16_t *indices, int numPrims);

extern Mat4 Model;
extern Mat4 View;
extern Mat4 Proj;
extern Mat4 ModelView;
extern Mat4 MVP;
extern int clipping;
