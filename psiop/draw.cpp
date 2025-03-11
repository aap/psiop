#include <Arduino.h>
#include <stdint.h>
#include <string.h>
#include "raster.h"

#include <math.h>
#include "vec4.h"

#include "draw.h"


/*
 * TODO:
 *	support texture coords??
 */

static rs::Canvas *canvas;

Mat4 Model;
Mat4 View;
Mat4 Proj;
Mat4 ModelView;
Mat4 MVP;
int clipping = 1;

Vec4 clipPlanes[6] = {
	{ 0.0f,  0.0f,  1.0f,  1.0f },
	{ 0.0f,  0.0f, -1.0f,  1.0f },
	{ 1.0f,  0.0f,  0.0f,  1.0f },
	{-1.0f,  0.0f,  0.0f,  1.0f },
	{ 0.0f,  1.0f,  0.0f,  1.0f },
	{ 0.0f, -1.0f,  0.0f,  1.0f }
};

static void drawPointCulled(VertCol *vin);
static void drawLineClipped(VertCol *vin);
static void drawTriClipped(VertCol *vin);

void
initDraw(uint16_t *fb, uint32_t *zb, int w, int h)
{
	static rs::Canvas CANV;
	canvas = &CANV;
	canvas->fb = (rs::u8*)fb;
	canvas->zbuf = zb;
	canvas->w = w;
	canvas->h = h;
	canvas->d = 0xFFFFFF - 1;

	clear(1, 0);
}

void
clear(int zbuf, uint16_t col)
{
	using namespace rs;

	int n = canvas->w*canvas->h;
	u16 *fb = (u16*)canvas->fb;
	while(n--) *fb++ = col;
	if(zbuf)
		memset(canvas->zbuf, 0, canvas->w*canvas->h*4);

	srScissorX0 = 0;
	srScissorX1 = canvas->w-1;
	srScissorY0 = 0;
	srScissorY1 = canvas->h-1;
}


void
updatemat(void)
{
	ModelView = m4mul(View, Model);
	MVP = m4mul(Proj, ModelView);
}

static Vec4 toscreen(Vec4 p) {
	p.x = (p.x+1.0f)/2.0f*canvas->w;
	p.y = (-p.y+1.0f)/2.0f*canvas->h;
	p.z = (-p.z+1.0f)/2.0f*canvas->d;
	return p;
}

static VertCol *storedVerts;

// we remember 1/z in w
static Vec4
project(Vec4 v)
{
	float q = 1.0f/v.w;
	v = toscreen(v4scale(q, v));
	v.w = q;
	return v;
}

static Vec4
xformProj(Vec4 v)
{
	return project(v4xform(&MVP, v));
}

static rs::Vertex
convVert(VertCol *vc)
{
	rs::Vertex v;
	float q = vc->pos.w;
	v.x = vc->pos.x * 16.0f;
	v.y = vc->pos.y * 16.0f;
	v.z = vc->pos.z;
	v.q = q;
	v.r = vc->r;
	v.g = vc->g;
	v.b = vc->b;
	v.a = vc->a;
	v.f = 0;
	v.s = 0.0f; //vc->u*q;
	v.t = 0.0f; //vc->v*q;
	return v;
}

void
drawPoint(VertCol v)
{
	if(clipping) {
		v.pos = v4xform(&MVP, v.pos);
		drawPointCulled(&v);
	} else {
		v.pos = xformProj(v.pos);
		rs::rasterPoint(canvas, convVert(&v));
	}
}

void
drawLine(VertCol v1, VertCol v2)
{
	if(clipping) {
		VertCol v[2];
		v[0] = v1;
		v[1] = v2;
		v[0].pos = v4xform(&MVP, v[0].pos);
		v[1].pos = v4xform(&MVP, v[1].pos);
		drawLineClipped(v);
	} else {
		v1.pos = xformProj(v1.pos);
		v2.pos = xformProj(v2.pos);
		rs::rasterLine(canvas, convVert(&v1), convVert(&v2));
	}
}

void
drawTriangle(VertCol v1, VertCol v2, VertCol v3)
{
	if(clipping) {
		VertCol v[3];
		v[0] = v1;
		v[1] = v2;
		v[2] = v3;
		v[0].pos = v4xform(&MVP, v[0].pos);
		v[1].pos = v4xform(&MVP, v[1].pos);
		v[2].pos = v4xform(&MVP, v[2].pos);
		drawTriClipped(v);
	} else {
		v1.pos = xformProj(v1.pos);
		v2.pos = xformProj(v2.pos);
		v3.pos = xformProj(v3.pos);
		rs::rasterTriangle(canvas, convVert(&v1), convVert(&v2), convVert(&v3));
	}
}

void
setVertices(VertCol *verts, int numVerts)
{
	storedVerts = verts;
	if(clipping) {
		for(int i = 0; i < numVerts; i++)
			verts[i].pos = v4xform(&MVP, verts[i].pos);
	} else {
		for(int i = 0; i < numVerts; i++)
			verts[i].pos = xformProj(verts[i].pos);
	}
}

void
drawIndexed(int prim, uint16_t *indices, int numPrims)
{
	VertCol v[3];

	switch(prim) {
	case 1:
		if(clipping) {
			while(numPrims--) {
				v[0] = storedVerts[indices[0]];
				drawPointCulled(v);

				indices += 1;
			}
		} else {
			while(numPrims--) {
				v[0] = storedVerts[indices[0]];
				rs::rasterPoint(canvas, convVert(&v[0]));

				indices += 1;
			}
		}
		break;

	case 2:
		if(clipping) {
			while(numPrims--) {
				v[0] = storedVerts[indices[0]];
				v[1] = storedVerts[indices[1]];
				drawLineClipped(v);

				indices += 2;
			}
		} else {
			while(numPrims--) {
				v[0] = storedVerts[indices[0]];
				v[1] = storedVerts[indices[1]];
				rs::rasterLine(canvas, convVert(&v[0]), convVert(&v[1]));

				indices += 2;
			}
		}
		break;

	case 3:
		if(clipping) {
			while(numPrims--) {
				v[0] = storedVerts[indices[0]];
				v[1] = storedVerts[indices[1]];
				v[2] = storedVerts[indices[2]];
				drawTriClipped(v);

				indices += 3;
			}
		} else {
			while(numPrims--) {
				v[0] = storedVerts[indices[0]];
				v[1] = storedVerts[indices[1]];
				v[2] = storedVerts[indices[2]];
				rs::rasterTriangle(canvas, convVert(&v[0]), convVert(&v[1]), convVert(&v[2]));

				indices += 3;
			}
		}
		break;
	}
}



/*
 *
 * Clipping
 *
 */



struct InterVertex
{
	Vec4 pos;
	Vec4 bary;
};

static void
interpVerts(rs::Vertex *out, InterVertex *clip, int n, VertCol *in)
{
	VertCol v;
	while(n--) {
		v.pos = clip->pos;

		v.r =	clip->bary.x*in[0].r +
			clip->bary.y*in[1].r +
			clip->bary.z*in[2].r;
		v.g =	clip->bary.x*in[0].g +
			clip->bary.y*in[1].g +
			clip->bary.z*in[2].g;
		v.b =	clip->bary.x*in[0].b +
			clip->bary.y*in[1].b +
			clip->bary.z*in[2].b;
		v.a =	clip->bary.x*in[0].a +
			clip->bary.y*in[1].a +
			clip->bary.z*in[2].a;

		v.pos = project(v.pos);
		*out = convVert(&v);

		clip++;
		out++;
	}
}

// sutherland-hodgman clipping
static int
clipPoly(InterVertex *in, int nin, InterVertex *out, Vec4 plane)
{
	int j;
	int nout;
	int x1, x2;
	float d1, d2;

	nout = 0;
	x1 = nin-1;
	for(j = 0; j < nin; j++) {
		x2 = j;
		d1 = v4dot(plane, in[x1].pos);
		d2 = v4dot(plane, in[x2].pos);
		if(d1 >= 0.0f)
			out[nout++] = in[x1];
		if(d1*d2 < 0.0f) {
			float inv = 1.0f/(d1 - d2);
			out[nout].pos = v4scale(inv, v4sub(v4scale(d1, in[x2].pos), v4scale(d2, in[x1].pos)));
			out[nout].bary = v4scale(inv, v4sub(v4scale(d1, in[x2].bary), v4scale(d2, in[x1].bary)));
			nout++;
		}
		x1 = x2;
	}
	return nout;
}

static int
clipLine(InterVertex *in, InterVertex *out, Vec4 plane)
{
	int nout;
	float d1, d2;

	nout = 0;
	d1 = v4dot(plane, in[0].pos);
	d2 = v4dot(plane, in[1].pos);
	if(d1 >= 0.0f)
		out[nout++] = in[0];
	if(d1*d2 < 0.0f) {
		float inv = 1.0f/(d1 - d2);
		out[nout].pos = v4scale(inv, v4sub(v4scale(d1, in[1].pos), v4scale(d2, in[0].pos)));
		out[nout].bary = v4scale(inv, v4sub(v4scale(d1, in[1].bary), v4scale(d2, in[0].bary)));
		nout++;
	}
	if(d2 >= 0.0f)
		out[nout++] = in[1];
	return nout != 2;
}

static void
drawPointCulled(VertCol *vin)
{
	Vec4 p = vin[0].pos;
	for(int i = 0; i < 6; i++)
		if(v4dot(clipPlanes[i], p) < 0.0f)
			return;
	vin[0].pos = project(vin[0].pos);
	rs::rasterPoint(canvas, convVert(&vin[0]));
}

static void
drawLineClipped(VertCol *vin)
{
	InterVertex buf[4];
	InterVertex *in = &buf[0];
	InterVertex *out = &buf[2];

	in[0].pos = vin[0].pos;
	in[1].pos = vin[1].pos;
	in[0].bary = vec4(1.0f, 0.0f, 0.0f, 0.0f);
	in[1].bary = vec4(0.0f, 1.0f, 0.0f, 0.0f);

	if(clipLine(in,  out, clipPlanes[0])) return;
	if(clipLine(out,  in, clipPlanes[1])) return;
	if(clipLine(in,  out, clipPlanes[2])) return;
	if(clipLine(out,  in, clipPlanes[3])) return;
	if(clipLine(in,  out, clipPlanes[4])) return;
	if(clipLine(out,  in, clipPlanes[5])) return;
	out = in;

	rs::Vertex verts[2];
	interpVerts(verts, out, 2, vin);
	rs::rasterLine(canvas, verts[0], verts[1]);
}

static void
drawTriClipped(VertCol *vin)
{
	InterVertex buf[18];
	InterVertex *in = &buf[0];
	InterVertex *out = &buf[9];

	in[0].pos = vin[0].pos;
	in[1].pos = vin[1].pos;
	in[2].pos = vin[2].pos;
	in[0].bary = vec4(1.0f, 0.0f, 0.0f, 0.0f);
	in[1].bary = vec4(0.0f, 1.0f, 0.0f, 0.0f);
	in[2].bary = vec4(0.0f, 0.0f, 1.0f, 0.0f);

	int nout = 3;
	if(nout = clipPoly(in,  nout, out, clipPlanes[0]), nout == 0) return;
	if(nout = clipPoly(out, nout,  in, clipPlanes[1]), nout == 0) return;
	if(nout = clipPoly(in,  nout, out, clipPlanes[2]), nout == 0) return;
	if(nout = clipPoly(out, nout,  in, clipPlanes[3]), nout == 0) return;
	if(nout = clipPoly(in,  nout, out, clipPlanes[4]), nout == 0) return;
	if(nout = clipPoly(out, nout,  in, clipPlanes[5]), nout == 0) return;
	out = in;

	rs::Vertex verts[9];
	interpVerts(verts, out, nout, vin);
	for(int i = 0; i < nout-2; i++)
		rs::rasterTriangle(canvas, verts[0], verts[i+1], verts[i+2]);
}
