namespace rs {

typedef int8_t i8;
typedef uint8_t u8;
typedef int16_t i16;
typedef uint16_t u16;
typedef int32_t i32;
typedef uint32_t u32;
typedef int64_t i64;
typedef uint64_t u64;



typedef struct Canvas Canvas;
struct Canvas
{
	u8 *fb;
	u32 *zbuf;
	int w, h, d;
};

typedef struct Texture Texture;
struct Texture          
{                       
	u8 *pixels;
	int w, h;
	int wrap;
};              

typedef struct Point3 Point3;
struct Point3
{
	int x, y, z;
};

typedef struct Color Color;
struct Color
{
	u8 r, g, b, a;
};

typedef struct Vertex Vertex;
struct Vertex
{
	i32 x, y, z;
	float q;	// 1/z
	u8 r, g, b, a;
	u8 f;		// fog
	float s, t;
};



Canvas *makecanvas(int w, int h);
Texture *maketexture(int w, int h);
void putpixel(Canvas *canvas, Point3 p, Color c);
void clearcanvas(Canvas *canvas);
void rasterPoint(Canvas *canvas, Vertex p);
void rasterLine(Canvas *canvas, Vertex p1, Vertex p2);
void rasterTriangle(Canvas *canvas, Vertex p1, Vertex p2, Vertex p3);

// not good
//void drawRect(Canvas *canvas, Point3 p1, Point3 p2, Color c);
//void drawLine(Canvas *canvas, Point3 p1, Point3 p2, Color c);

/*
 * Render States
 */
enum TextureWrap {
	WRAP_REPEAT,
	WRAP_CLAMP,
	WRAP_BORDER,
};

enum TextureFunction {
	TFUNC_MODULATE,
	TFUNC_DECAL,
	TFUNC_HIGHLIGHT,
	TFUNC_HIGHLIGHT2,
};

enum AlphaTestFunc {
	ALPHATEST_NEVER,
	ALPHATEST_ALWAYS,
	ALPHATEST_LESS,
	ALPHATEST_LEQUAL,
	ALPHATEST_EQUAL,
	ALPHATEST_GEQUAL,
	ALPHATEST_GREATER,
	ALPHATEST_NOTEQUAL,
};

enum AlphaTestFail {
	ALPHAFAIL_KEEP,
	ALPHAFAIL_FB_ONLY,
	ALPHAFAIL_ZB_ONLY,
};

enum DepthTestFunc {
	DEPTHTEST_NEVER,
	DEPTHTEST_ALWAYS,
	DEPTHTEST_GEQUAL,
	DEPTHTEST_GREATER,
};

// The blend equation is
// out = ((A - B) * C >> 7) + D
// A, B and D select the color, C the alpha value
enum AlphaBlendOp {
	ALPHABLEND_SRC,
	ALPHABLEND_DST,
	ALPHABLEND_ZERO,
	ALPHABLEND_FIX = ALPHABLEND_ZERO,
};

extern int srScissorX0, srScissorX1;
extern int srScissorY0, srScissorY1;
extern int srDepthTestEnable;
extern int srDepthTestFunction;
extern int srWriteZ;
extern int srAlphaTestEnable;
extern int srAlphaTestFunction;
extern int srAlphaTestReference;
extern int srAlphaTestFail;
extern int srAlphaBlendEnable;
extern int srAlphaBlendA;
extern int srAlphaBlendB;
extern int srAlphaBlendC;
extern int srAlphaBlendD;
extern int srAlphaBlendFix;
extern int srTexEnable;
extern Texture *srTexture;
extern int srWrapU;
extern int srWrapV;
extern Color srBorder;
extern int srTexUseAlpha;
extern int srTexFunc;
extern int srFogEnable;
extern Color srFogCol;

}
