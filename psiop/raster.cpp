#include <stdint.h>

// should get rid of malloc
#include <stdlib.h>
#include <string.h>

//#define assert(X)

#include "raster.h"

namespace rs {

#define CEIL(p) (((p)+15) >> 4)

// render states
int srScissorX0, srScissorX1;
int srScissorY0, srScissorY1;
int srDepthTestEnable = 1;
int srDepthTestFunction = DEPTHTEST_GEQUAL;
int srWriteZ = 1;
int srAlphaTestEnable = 1;
int srAlphaTestFunction = ALPHATEST_ALWAYS;
int srAlphaTestReference;
int srAlphaTestFail = ALPHAFAIL_FB_ONLY;
int srAlphaBlendEnable = 1;
int srAlphaBlendA = ALPHABLEND_SRC;
int srAlphaBlendB = ALPHABLEND_DST;
int srAlphaBlendC = ALPHABLEND_SRC;
int srAlphaBlendD = ALPHABLEND_DST;
int srAlphaBlendFix = 0x80;
int srTexEnable = 0;
Texture *srTexture;
int srWrapU = WRAP_REPEAT;
int srWrapV = WRAP_REPEAT;
Color srBorder = { 255, 0, 0, 255 };
int srTexUseAlpha = 1;
int srTexFunc = TFUNC_MODULATE;
int srFogEnable = 0;
Color srFogCol = { 0, 0, 0, 0 };

int clamp(int x) { if(x < 0) return 0; if(x > 255) return 255; return x; }
int abs(int x) { return x < 0 ? -x : x; }

#if 0
Canvas*
makecanvas(int w, int h)
{
	Canvas *canv;
	canv = (Canvas*)malloc(sizeof(*canv) + w*h*(4+4));
if(canv == NULL) {
	for(;;);
}
	canv->w = w;
	canv->h = h;
	canv->fb = ((u8*)canv + sizeof(*canv));
	canv->zbuf = (u32*)(canv->fb + w*h*4);
	return canv;
}

Texture*
maketexture(int w, int h)
{
	Texture *t;
	t = (Texture*)malloc(sizeof(*t) + w*h*4);
	t->w = w;
	t->h = h;
	t->pixels = (u8*)t + sizeof(*t);
	t->wrap = 0x11; // wrap u and v
	return t;
}

void
clearcanvas(Canvas *canvas)
{
#if 0
	memset(canvas->fb, 0, canvas->w*canvas->h*4);
#else
	memset(canvas->fb, 0, canvas->w*canvas->h*2);
#endif
	memset(canvas->zbuf, 0, canvas->w*canvas->h*4);
}
#endif

void
writefb(Canvas *canvas, int x, int y, Color c)
{
#if 0
	u8 *px = &canvas->fb[(y*canvas->w + x)*4];
	u32 *z = &canvas->zbuf[y*canvas->w + x];

	px[3] = c.r;
	px[2] = c.g;
	px[1] = c.b;
	px[0] = c.a;
#else
	u16 *fb = (u16*)canvas->fb;
	int r = c.r >> 3;
	int g = c.g >> 2;
	int b = c.b >> 3;
	fb[y*canvas->w + x] = (r<<11) | (g<<5) | b;
#endif
}

void
putpixel(Canvas *canvas, Point3 p, Color c)
{
	// scissor test
	if(p.x < srScissorX0 || p.x > srScissorX1 ||
	   p.y < srScissorY0 || p.y > srScissorY1)
		return;

//	u8 *px = &canvas->fb[(p.y*canvas->w + p.x)*4];
	u32 *z = &canvas->zbuf[p.y*canvas->w + p.x];

	int fbwrite = 1;
	int zbwrite = srWriteZ;

	// alpha test
	if(srAlphaTestEnable){
		int fail;
		switch(srAlphaTestFunction){
		case ALPHATEST_NEVER:
			fail = 1;
			break;
		case ALPHATEST_ALWAYS:
		default:
			fail = 0;
			break;
		case ALPHATEST_LESS:
			fail = c.a >= srAlphaTestReference;
			break;
		case ALPHATEST_LEQUAL:
			fail = c.a > srAlphaTestReference;
			break;
		case ALPHATEST_EQUAL:
			fail = c.a != srAlphaTestReference;
			break;
		case ALPHATEST_GEQUAL:
			fail = c.a < srAlphaTestReference;
			break;
		case ALPHATEST_GREATER:
			fail = c.a <= srAlphaTestReference;
			break;
		case ALPHATEST_NOTEQUAL:
			fail = c.a == srAlphaTestReference;
			break;
		}
		if(fail){
			switch(srAlphaTestFail){
			case ALPHAFAIL_KEEP:
				return;
			case ALPHAFAIL_FB_ONLY:
				zbwrite = 0;
				break;
			case ALPHAFAIL_ZB_ONLY:
				fbwrite = 0;
			}
		}
	}

	// ztest
	if(srDepthTestEnable){
		switch(srDepthTestFunction){
		case DEPTHTEST_NEVER:
			return;
		case DEPTHTEST_ALWAYS:
			break;
		case DEPTHTEST_GEQUAL:
			if((u32)p.z < *z)
				return;
			break;
		case DEPTHTEST_GREATER:
			if((u32)p.z <= *z)
				return;
			break;
		}
	}

// no blending for now
#if 0
	Color d = { px[3], px[2], px[1], px[0] };

	// blend
	if(srAlphaBlendEnable){
		int ar, ag, ab;
		int br, bg, bb;
		int dr, dg, db;
		int ca;
		switch(srAlphaBlendA){
		case ALPHABLEND_SRC:
			ar = c.r;
			ag = c.g;
			ab = c.b;
			break;
		case ALPHABLEND_DST:
			ar = d.r;
			ag = d.g;
			ab = d.b;
			break;
		case ALPHABLEND_ZERO:
			ar = 0;
			ag = 0;
			ab = 0;
			break;
		default: assert(0);
		}
		switch(srAlphaBlendB){
		case ALPHABLEND_SRC:
			br = c.r;
			bg = c.g;
			bb = c.b;
			break;
		case ALPHABLEND_DST:
			br = d.r;
			bg = d.g;
			bb = d.b;
			break;
		case ALPHABLEND_ZERO:
			br = 0;
			bg = 0;
			bb = 0;
			break;
		default: assert(0);
		}
		switch(srAlphaBlendC){
		case ALPHABLEND_SRC:
			ca = c.a;
			break;
		case ALPHABLEND_DST:
			ca = d.a;
			break;
		case ALPHABLEND_FIX:
			ca = srAlphaBlendFix;
			break;
		default: assert(0);
		}
		switch(srAlphaBlendD){
		case ALPHABLEND_SRC:
			dr = c.r;
			dg = c.g;
			db = c.b;
			break;
		case ALPHABLEND_DST:
			dr = d.r;
			dg = d.g;
			db = d.b;
			break;
		case ALPHABLEND_ZERO:
			dr = 0;
			dg = 0;
			db = 0;
			break;
		default: assert(0);
		}

		int r, g, b;
		r = ((ar - br) * ca >> 7) + dr;
		g = ((ag - bg) * ca >> 7) + dg;
		b = ((ab - bb) * ca >> 7) + db;

		c.r = clamp(r);
		c.g = clamp(g);
		c.b = clamp(b);
	}
#endif

	if(fbwrite)
		writefb(canvas, p.x, p.y, c);
	if(zbwrite)
		*z = p.z;
}

Color
sampletex_nearest(int u, int v)
{
	Texture *tex = srTexture;

	const int usize = tex->w;
	const int vsize = tex->h;

	int iu = u >> 4;
	int iv = v >> 4;

	switch(srWrapU){
	case WRAP_REPEAT:
		iu %= usize;
		break;
	case WRAP_CLAMP:
		if(iu < 0) iu = 0;
		if(iu >= usize) iu = usize-1;
		break;
	case WRAP_BORDER:
		if(iu < 0 || iu >= usize)
			return srBorder;
	}

	switch(srWrapV){
	case WRAP_REPEAT:
		iv %= vsize;
		break;
	case WRAP_CLAMP:
		if(iv < 0) iv = 0;
		if(iv >= vsize) iv = vsize-1;
		break;
	case WRAP_BORDER:
		if(iv < 0 || iv >= vsize)
			return srBorder;
	}

	u8 *cp = &tex->pixels[(iv*tex->w + iu)*4];
	Color c = { cp[0], cp[1], cp[2], cp[3] };
	return c;
}

// t is texture, f is fragment
Color
texfunc(Color t, Color f)
{
	int r, g, b, a;
	switch(srTexFunc){
	case TFUNC_MODULATE:
		r = t.r * f.r >> 7;
		g = t.g * f.g >> 7;
		b = t.b * f.b >> 7;
		a = srTexUseAlpha ?
			t.a * f.a >> 7 :
			f.a;
		break;
	case TFUNC_DECAL:
	default:
		r = t.r;
		g = t.g;
		b = t.b;
		a = srTexUseAlpha ? t.a : f.a;
		break;
	case TFUNC_HIGHLIGHT:
		r = (t.r * f.r >> 7) + f.a;
		g = (t.g * f.g >> 7) + f.a;
		b = (t.b * f.b >> 7) + f.a;
		a = srTexUseAlpha ?
			t.a + f.a :
			f.a;
		break;
	case TFUNC_HIGHLIGHT2:
		r = (t.r * f.r >> 7) + f.a;
		g = (t.g * f.g >> 7) + f.a;
		b = (t.b * f.b >> 7) + f.a;
		a = srTexUseAlpha ? t.a : f.a;
		break;
	}
	Color v;
	v.r = clamp(r);
	v.g = clamp(g);
	v.b = clamp(b);
	v.a = clamp(a);
	return v;
}

Point3 mkpnt(int x, int y, int z) { Point3 p = { x, y, z}; return p; }

#if 0
void
drawRect(Canvas *canvas, Point3 p1, Point3 p2, Color c)
{
	int x, y;
	for(y = p1.y; y <= p2.y; y++)
		for(x = p1.x; x <= p2.x; x++)
			putpixel(canvas, mkpnt(x, y, 0), c);
}

void
drawLine(Canvas *canvas, Point3 p1, Point3 p2, Color c)
{
	int dx, dy;
	int incx, incy;
	int e;
	int x, y;

	dx = abs(p2.x-p1.x);
	incx = p2.x > p1.x ? 1 : -1;
	dy = abs(p2.y-p1.y);
	incy = p2.y > p1.y ? 1 : -1;
	e = 0;
	if(dx == 0){
		for(y = p1.y; y != p2.y; y += incy)
			putpixel(canvas, mkpnt(p1.x, y, 0), c);
	}else if(dx > dy){
		y = p1.y;
		for(x = p1.x; x != p2.x; x += incx){
			putpixel(canvas, mkpnt(x, y, 0), c);
			e += dy;
			if(2*e >= dx){
				e -= dx;
				y += incy;
			}
		}
	}else{
		x = p1.x;
		for(y = p1.y; y != p2.y; y += incy){
			putpixel(canvas, mkpnt(x, y, 0), c);
			e += dx;
			if(2*e >= dy){
				e -= dy;
				x += incx;
			}
		}
	}
}
#endif

Color
shadePixel(Color c, float s, float t, float q, int f)
{
	if(srTexEnable && srTexture){
		float w = 1.0f/q;
		s = s * w;
		t = t * w;
		int u = s * srTexture->w * 16;
		int v = t * srTexture->h * 16;
		Color texc = sampletex_nearest(u, v);
		c = texfunc(texc, c);
	}
	if(srFogEnable){
		c.r = (f*c.r >> 8) + ((255 - f)*srFogCol.r >> 8);
		c.g = (f*c.g >> 8) + ((255 - f)*srFogCol.g >> 8);
		c.b = (f*c.b >> 8) + ((255 - f)*srFogCol.b >> 8);
	}
	return c;
}

void
rasterPoint(Canvas *canvas, Vertex v)
{
	Color c = { v.r, v.g, v.b, v.a };
	putpixel(canvas, mkpnt((v.x+8)>>4, (v.y+8)>>4, v.z),
		shadePixel(c, v.s, v.t, v.q, v.f));
}

/*
	attibutes we want to interpolate:
	R G B A
	U V / S T Q
	X Y Z F
*/

struct PrimAttribs
{
	i64 z;
	i32 r, g, b, a;
	i32 f;
	float s, t;
	float q;
};

static void
add1(struct PrimAttribs *a, struct PrimAttribs *b)
{
	a->z += b->z;
	a->r += b->r;
	a->g += b->g;
	a->b += b->b;
	a->a += b->a;
	a->f += b->f;
	a->s += b->s;
	a->t += b->t;
	a->q += b->q;
}

static void
sub1(struct PrimAttribs *a, struct PrimAttribs *b)
{
	a->z -= b->z;
	a->r -= b->r;
	a->g -= b->g;
	a->b -= b->b;
	a->a -= b->a;
	a->f -= b->f;
	a->s -= b->s;
	a->t -= b->t;
	a->q -= b->q;
}

static void
guard(struct PrimAttribs *a)
{
	if(a->z < 0) a->z = 0;
	else if(a->z > 0x3FFFFFFFC000LL) a->z = 0x3FFFFFFFC000LL;
	if(a->r < 0) a->r = 0;
	else if(a->r > 0xFF000) a->r = 0xFF000;
	if(a->g < 0) a->g = 0;
	else if(a->g > 0xFF000) a->g = 0xFF000;
	if(a->b < 0) a->b = 0;
	else if(a->b > 0xFF000) a->b = 0xFF000;
	if(a->a < 0) a->a = 0;
	else if(a->a > 0xFF000) a->a = 0xFF000;
	if(a->f < 0) a->f = 0;
	else if(a->f > 0xFF000) a->f = 0xFF000;
}

struct RasTri
{
	int x, y;
	int ymid, yend;
	int right;
	int e[2], dx[3], dy[3];
	struct PrimAttribs gx, gy, v, s;
};

static int
triangleSetup(struct RasTri *tri, Vertex v1, Vertex v2, Vertex v3)
{
	int dx1, dx2, dx3;
	int dy1, dy2, dy3;

	dy1 = v3.y - v1.y;	// long edge
	if(dy1 == 0) return 1;
	dx1 = v3.x - v1.x;
	dx2 = v2.x - v1.x;	// first small edge
	dy2 = v2.y - v1.y;
	dx3 = v3.x - v2.x;	// second small edge
	dy3 = v3.y - v2.y;

	// this is twice the triangle area
	const int area = dx2*dy1 - dx1*dy2;
	if(area == 0) return 1;
	// figure out if 0 or 1 is the right edge
	tri->right = area < 0;

	/* The gradients are to step whole pixels,
	 * so they are pre-multiplied by 16. */

	float denom = 16.0f/area;
	// gradients x
#define GX(p) ((v2.p - v1.p)*dy1 - (v3.p - v1.p)*dy2)
	tri->gx.z = (((i64)v2.z - (i64)v1.z)*dy1 - ((i64)v3.z - (i64)v1.z)*dy2)*denom * 16384;
	tri->gx.r = GX(r)*denom * 4096;
	tri->gx.g = GX(g)*denom * 4096;
	tri->gx.b = GX(b)*denom * 4096;
	tri->gx.a = GX(a)*denom * 4096;
	tri->gx.f = GX(f)*denom * 4096;
	tri->gx.s = GX(s)*denom;
	tri->gx.t = GX(t)*denom;
	tri->gx.q = GX(q)*denom;

	// gradients y
	denom = -denom;
#define GY(p) ((v2.p - v1.p)*dx1 - (v3.p - v1.p)*dx2)
	tri->gy.z = (((i64)v2.z - (i64)v1.z)*dx1 - ((i64)v3.z - (i64)v1.z)*dx2)*denom * 16384;
	tri->gy.r = GY(r)*denom * 4096;
	tri->gy.g = GY(g)*denom * 4096;
	tri->gy.b = GY(b)*denom * 4096;
	tri->gy.a = GY(a)*denom * 4096;
	tri->gy.f = GY(f)*denom * 4096;
	tri->gy.s = GY(s)*denom;
	tri->gy.t = GY(t)*denom;
	tri->gy.q = GY(q)*denom;

	tri->ymid = CEIL(v2.y);
	tri->yend = CEIL(v3.y);

	tri->y = CEIL(v1.y);
	tri->x = CEIL(v1.x);

	tri->dy[0] = dy2<<4;	// upper edge
	tri->dy[1] = dy1<<4;	// lower edge
	tri->dy[2] = dy3<<4;	// long edge
	tri->dx[0] = dx2<<4;
	tri->dx[1] = dx1<<4;
	tri->dx[2] = dx3<<4;

	// prestep to land on pixel center

	int stepx = v1.x - (tri->x<<4);
	int stepy = v1.y - (tri->y<<4);
	tri->e[0] = (-stepy*tri->dx[0] + stepx*tri->dy[0]) >> 4;
	tri->e[1] = (-stepy*tri->dx[1] + stepx*tri->dy[1]) >> 4;

	// attributes along interpolated edge	
	tri->v.z = (i64)v1.z*16384 - (stepy*tri->gy.z + stepx*tri->gx.z)/16;
	tri->v.r = v1.r*4096 - (stepy*tri->gy.r + stepx*tri->gx.r)/16;
	tri->v.g = v1.g*4096 - (stepy*tri->gy.g + stepx*tri->gx.g)/16;
	tri->v.b = v1.b*4096 - (stepy*tri->gy.b + stepx*tri->gx.b)/16;
	tri->v.a = v1.a*4096 - (stepy*tri->gy.a + stepx*tri->gx.a)/16;
	tri->v.f = v1.f*4096 - (stepy*tri->gy.f + stepx*tri->gx.f)/16;
	tri->v.s = v1.s - (stepy*tri->gy.s + stepx*tri->gx.s)/16.0f;
	tri->v.t = v1.t - (stepy*tri->gy.t + stepx*tri->gx.t)/16.0f;
	tri->v.q = v1.q - (stepy*tri->gy.q + stepx*tri->gx.q)/16.0f;

	return 0;
}

void
rasterTriangle(Canvas *canvas, Vertex v1, Vertex v2, Vertex v3)
{
	Color c;
	struct RasTri tri;
	int stepx, stepy;

	// Sort such that we have from top to bottom v1,v2,v3
	if(v2.y < v1.y){ Vertex tmp = v1; v1 = v2; v2 = tmp; }
	if(v3.y < v1.y){ Vertex tmp = v1; v1 = v3; v3 = tmp; }
	if(v3.y < v2.y){ Vertex tmp = v2; v2 = v3; v3 = tmp; }

	if(triangleSetup(&tri, v1, v2, v3))
		return;

	// Current scanline start and end
	int xn[2] = { tri.x, tri.x };
	int a = !tri.right;	// left edge
	int b = tri.right;	// right edge

	// If upper triangle has no height, only do the lower part
	if(tri.dy[0] == 0)
		goto secondtri;
	while(tri.y < tri.yend){
		/* TODO: is this the right way to step the edges? */

		/* Step x and interpolated value down left edge */
		while(tri.e[a] <= -tri.dy[a]){
			xn[a]--;
			tri.e[a] += tri.dy[a];
			sub1(&tri.v, &tri.gx);
		}
		while(tri.e[a] > 0){
			xn[a]++;
			tri.e[a] -= tri.dy[a];
			add1(&tri.v, &tri.gx);
		}

		/* Step x down right edge */
		while(tri.e[b] <= -tri.dy[b]){
			xn[b]--;
			tri.e[b] += tri.dy[b];
		}
		while(tri.e[b] > 0){
			xn[b]++;
			tri.e[b] -= tri.dy[b];
		}

		// When we reach the mid vertex, change state and jump to start of loop again
		// TODO: this is a bit ugly in here...can we fix it?
		if(tri.y == tri.ymid){
		secondtri:
			tri.dx[0] = tri.dx[2];
			tri.dy[0] = tri.dy[2];
			// Either the while prevents this or we returned early because dy1 == 0
			assert(tri.dy[0] != 0);
			stepx = v2.x - (xn[0]<<4);
			stepy = v2.y - (tri.y<<4);
			tri.e[0] = (-stepy*tri.dx[0] + stepx*tri.dy[0]) >> 4;

			tri.ymid = -1;	// so we don't do this again
			continue;
		}

		/* Rasterize one line */
		tri.s = tri.v;
		for(tri.x = xn[a]; tri.x < xn[b]; tri.x++){
			guard(&tri.s);
			c.r = tri.s.r >> 12;
			c.g = tri.s.g >> 12;
			c.b = tri.s.b >> 12;
			c.a = tri.s.a >> 12;
			putpixel(canvas, mkpnt(tri.x, tri.y, tri.s.z>>14),
				shadePixel(c, tri.s.s, tri.s.t, tri.s.q, tri.s.f>>12));
			add1(&tri.s, &tri.gx);
		}

		/* Step in y */
		tri.y++;
		tri.e[a] += tri.dx[a];
		tri.e[b] += tri.dx[b];
		add1(&tri.v, &tri.gy);
	}
}



void
rasterLine(Canvas *canvas, Vertex v1, Vertex v2)
{
	int x = (v1.x+8)>>4;
	int y = (v1.y+8)>>4;
	int xe = (v2.x+8)>>4;
	int ye = (v2.y+8)>>4;
	int dx = v2.x-v1.x;
	int dy = v2.y-v1.y;
	int adx = abs(dx);
	int ady = abs(dy);
	int incx = dx > 0 ? 1 : -1;
	int incy = dy > 0 ? 1 : -1;
	int e;
	struct PrimAttribs gx, gy, v, s;

	int stepx = incx*(v1.x - (x<<4));
	int stepy = incy*(v1.y - (y<<4));

	float lsq = dx*dx + dy*dy;
	if(lsq == 0.0f)
		return;

	float scl = 16.0f*adx/lsq;
	gx.z = ((i64)v2.z - (i64)v1.z)*scl * 16384;
	gx.r = (v2.r - v1.r)*scl * 4096;
	gx.g = (v2.g - v1.g)*scl * 4096;
	gx.b = (v2.b - v1.b)*scl * 4096;
	gx.a = (v2.a - v1.a)*scl * 4096;
	gx.f = (v2.f - v1.f)*scl * 4096;
	gx.s = (v2.s - v1.s)*scl;
	gx.t = (v2.t - v1.t)*scl;
	gx.q = (v2.q - v1.q)*scl;

	scl = 16.0f*ady/lsq;
	gy.z = ((i64)v2.z - (i64)v1.z)*scl * 16384;
	gy.r = (v2.r - v1.r)*scl * 4096;
	gy.g = (v2.g - v1.g)*scl * 4096;
	gy.b = (v2.b - v1.b)*scl * 4096;
	gy.a = (v2.a - v1.a)*scl * 4096;
	gy.f = (v2.f - v1.f)*scl * 4096;
	gy.s = (v2.s - v1.s)*scl;
	gy.t = (v2.t - v1.t)*scl;
	gy.q = (v2.q - v1.q)*scl;

	v.z = (i64)v1.z*16384 - (stepy*gy.z + stepx*gx.z)/16;
	v.r = v1.r*4096 - (stepy*gy.r + stepx*gx.r)/16;
	v.g = v1.g*4096 - (stepy*gy.g + stepx*gx.g)/16;
	v.b = v1.b*4096 - (stepy*gy.b + stepx*gx.b)/16;
	v.a = v1.a*4096 - (stepy*gy.a + stepx*gx.a)/16;
	v.f = v1.f*4096 - (stepy*gy.f + stepx*gx.f)/16;
	v.s = v1.s - (stepy*gy.s + stepx*gx.s)/16.0f;
	v.t = v1.t - (stepy*gy.t + stepx*gx.t)/16.0f;
	v.q = v1.q - (stepy*gy.q + stepx*gx.q)/16.0f;

	Color c;
	if(dx == 0) {
		for(; y != ye; y += incy) {
			s = v;
			guard(&s);
			c = (Color){ s.r>>12, s.g>>12, s.b>>12, s.a>>12 };
			putpixel(canvas, (Point3){x, y, s.z>>14},
				shadePixel(c, s.s, s.t, s.q, s.f>>12));
			add1(&v, &gy);
		}
	} else if(dy == 0) {
		for(; x != xe; x += incx) {
			s = v;
			guard(&s);
			c = (Color){ s.r>>12, s.g>>12, s.b>>12, s.a>>12 };
			putpixel(canvas, (Point3){x, y, s.z>>14},
				shadePixel(c, s.s, s.t, s.q, s.f>>12));
			add1(&v, &gx);
		}
	} else if(adx < ady) {
		e = (-stepy*adx + stepx*ady) >> 4;
		for(; y != ye; y += incy) {
			if(2*e > ady) {
				e -= ady;
				x += incx;
				add1(&v, &gx);
			}
			s = v;
			guard(&s);
			c = (Color){ s.r>>12, s.g>>12, s.b>>12, s.a>>12 };
			putpixel(canvas, (Point3){x, y, s.z>>14},
				shadePixel(c, s.s, s.t, s.q, s.f>>12));
			e += adx;
			add1(&v, &gy);
		}
	} else {
		e = (stepy*adx - stepx*ady) >> 4;
		for(; x != xe; x += incx) {
			if(2*e > adx) {
				e -= adx;
				y += incy;
				add1(&v, &gy);
			}
			s = v;
			guard(&s);
			c = (Color){ s.r>>12, s.g>>12, s.b>>12, s.a>>12 };
			putpixel(canvas, (Point3){x, y, s.z>>14},
				shadePixel(c, s.s, s.t, s.q, s.f>>12));
			e += ady;
			add1(&v, &gx);
		}
	}
}

}
