#include <Adafruit_GFX.h>    // Core graphics library
#include <Adafruit_ST7735.h> // Hardware-specific library for ST7735
#include <Adafruit_ST7789.h> // Hardware-specific library for ST7789
#include <SPI.h>

//#include "cpx.h"
#include "quat.h"
#include "vec4.h"

#include "raster.h"

#define nelem(a) (sizeof(a)/sizeof(*a))
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t u8;

#include "draw.h"

//#define RES160


#define BLACK ST77XX_BLACK
#define WHITE ST77XX_WHITE
#define RED   ST77XX_RED
#define GREEN ST77XX_GREEN
#define BLUE  ST77XX_BLUE

#define TFT_CS       5
#define TFT_RST      4
#define TFT_DC       2
Adafruit_ST7735 tft = Adafruit_ST7735(TFT_CS, TFT_DC, TFT_RST);

#ifdef RES160
#define WIDTH 160
#define HEIGHT 128
#define DISPTYPE INITR_GREENTAB
#else
#define WIDTH 128
#define HEIGHT 128
#define DISPTYPE INITR_144GREENTAB
#endif

#define ASPECT ((float)WIDTH/HEIGHT)


#ifndef PI
#define PI 3.14159f
#endif
#define TAU (2.0f*(float)PI)


u16 *framebuffer;
u32 *zbuffer;

void
sendfb(void)
{
	tft.startWrite();
	tft.setAddrWindow(0, 0, WIDTH, HEIGHT);
	tft.writePixels(framebuffer, WIDTH*HEIGHT);
	tft.endWrite();
}


float map(float x, float x1,float x2, float y1,float y2) { return (x-x1)/(x2-x1) * (y2-y1) + y1; }


Mat4 quat2mat4(Quat r) {
	Quat rx = qsandwich(r, quat(0.0f, 1.0f, 0.0f, 0.0f));
	Quat ry = qsandwich(r, quat(0.0f, 0.0f, 1.0f, 0.0f));
	Quat rz = qsandwich(r, quat(0.0f, 0.0f, 0.0f, 1.0f));
	return (Mat4){
		{ rx.x, rx.y, rx.z, 0.0f },
		{ ry.x, ry.y, ry.z, 0.0f },
		{ rz.x, rz.y, rz.z, 0.0f },
		{ 0.0f, 0.0f, 0.0f, 1.0f }
	};
}




static struct {
	int close;
	int havepoint;
	VertCol firstp;
	VertCol lastp;
} linestate;

void
beginLine(int close)
{
	linestate.close = close;
	linestate.havepoint = 0;
}

void
endLine(void)
{
	if(linestate.havepoint && linestate.close)
		drawLine(linestate.firstp, linestate.lastp);
	linestate.havepoint = 0;
}

void
vertex(Vec4 p)
{
	VertCol v = { p, 255, 255, 255, 255 };
	if(linestate.havepoint) {
		drawLine(linestate.lastp, v);
		linestate.lastp = v;
	} else {
		linestate.havepoint = 1;
		linestate.firstp = v;
		linestate.lastp = v;
	}
}


void
cube(void)
{
	VertCol p[8] = {
		{ vec4( 1.0f,  1.0f,  1.0f, 1.0f),   0xFF, 0xFF, 0xFF, 0xFF },
		{ vec4( 1.0f,  1.0f, -1.0f, 1.0f),   0xFF, 0xFF, 0xFF, 0xFF },
		{ vec4( 1.0f, -1.0f,  1.0f, 1.0f),   0xFF, 0xFF, 0xFF, 0xFF },
		{ vec4( 1.0f, -1.0f, -1.0f, 1.0f),   0xFF, 0xFF, 0xFF, 0xFF },
		{ vec4(-1.0f,  1.0f,  1.0f, 1.0f),   0xFF, 0xFF, 0xFF, 0xFF },
		{ vec4(-1.0f,  1.0f, -1.0f, 1.0f),   0xFF, 0xFF, 0xFF, 0xFF },
		{ vec4(-1.0f, -1.0f,  1.0f, 1.0f),   0xFF, 0xFF, 0xFF, 0xFF },
		{ vec4(-1.0f, -1.0f, -1.0f, 1.0f),   0xFF, 0xFF, 0xFF, 0xFF },
	};
	static u16 lines[] = {
		0,1,
		1,3,
		3,2,
		2,0,

		4,5,
		5,7,
		7,6,
		6,4,

		0,4,
		1,5,
		3,7,
		2,6,
	};
	float s = 0.8f;
	Model = m4mul(m4scale(s, s, s), Model);
	updatemat();
	setVertices(p, nelem(p));
	drawIndexed(2, lines, nelem(lines)/2);
}

void
filled_cube(void)
{
	VertCol p[8] = {
		{ vec4( 1.0f,  1.0f,  1.0f, 1.0f),   0xFF, 0xFF, 0xFF, 0xFF },
		{ vec4( 1.0f,  1.0f, -1.0f, 1.0f),   0xFF, 0xFF, 0x00, 0xFF },
		{ vec4( 1.0f, -1.0f,  1.0f, 1.0f),   0xFF, 0x00, 0xFF, 0xFF },
		{ vec4( 1.0f, -1.0f, -1.0f, 1.0f),   0xFF, 0x00, 0x00, 0xFF },
		{ vec4(-1.0f,  1.0f,  1.0f, 1.0f),   0x00, 0xFF, 0xFF, 0xFF },
		{ vec4(-1.0f,  1.0f, -1.0f, 1.0f),   0x00, 0xFF, 0x00, 0xFF },
		{ vec4(-1.0f, -1.0f,  1.0f, 1.0f),   0x00, 0x00, 0xFF, 0xFF },
		{ vec4(-1.0f, -1.0f, -1.0f, 1.0f),   0x00, 0x00, 0x00, 0xFF },
	};
	static u16 tris[] = {
		// right
		0, 1, 2,
		1, 2, 3,

		// left
		4, 5, 6,
		5, 6, 7,

		// top
		0, 2, 4,
		2, 4, 6,

		// bottom
		1, 3, 5,
		3, 5, 7,

		// front
		0, 1, 4, 
		1, 4, 5,

		// back
		2, 3, 6, 
		3, 6, 7
	};
	float s = 0.8f;
	Model = m4mul(m4scale(s, s, s), Model);
	updatemat();
	setVertices(p, nelem(p));
	drawIndexed(3, tris, nelem(tris)/3);
}

void
axes(void)
{
	Model = m4ident();
	updatemat();
	VertCol p01 = { vec4(0.0f, 0.0f, 0.0f, 1.0f), 255, 0, 0, 255 };
	VertCol p1  = { vec4(1.0f, 0.0f, 0.0f, 1.0f), 255, 0, 0, 255 };
	VertCol p02 = { vec4(0.0f, 0.0f, 0.0f, 1.0f), 0, 255, 0, 255 };
	VertCol p2  = { vec4(0.0f, 1.0f, 0.0f, 1.0f), 0, 255, 0, 255 };
	VertCol p03 = { vec4(0.0f, 0.0f, 0.0f, 1.0f), 0, 0, 255, 255 };
	VertCol p3  = { vec4(0.0f, 0.0f, 1.0f, 1.0f), 0, 0, 255, 255 };
	drawLine(p01, p1);
	drawLine(p02, p2);
	drawLine(p03, p3);
}

Vec4 sph2cart(float r, float theta, float phi) {
	return vec4(r*cos(phi)*cos(theta),
	            r*sin(phi)*cos(theta),
	            r*sin(theta), 1.0f);
}

void
deathstar(void)
{
	float r = 0.85f;
	VertCol v, u;

	Model = m4mul(m4scale(1.0f, 1.0f, 0.85f), Model);
	updatemat();

	v.r = v.g = v.b = v.a = 255;
	u.r = u.g = u.b = u.a = 255;
	// sphere
	for(int j = 1; j < 22; j++) {
		float theta = map(j, 0,22, TAU/4,-TAU/2);
		for(int i = 0; i < 34; i++) {
			float phi = map(i, 0,34, 0,TAU);
			Vec4 p = sph2cart(r, theta, phi);
			float dsq = p.y*p.y + p.z*p.z;
			if(dsq > sq(r*0.35f) || p.x < 0) {
				v.pos = p;
				drawPoint(v);
			}
		}
	}

	// equator
	beginLine(0);
	for(int i = 2; i < 33; i++) {
		float phi = map(i, 0,34, 0,TAU);
		Vec4 p = sph2cart(r, 0.0f, phi);
		vertex(p);
	}
	endLine();

	// laser - outer ring
	beginLine(1);
	for(int i = 0; i < 15; i++) {
		float phi = map(i, 0,15, 0,TAU);
		Vec4 p = sph2cart(r, TAU/4-0.37f, phi);
		vertex(vec4(p.z, p.y, p.x, 1.0f));
	}
	endLine();

	// laser - inner ring
	beginLine(1);
	for(int i = 0; i < 15; i++) {
		float phi = map(i, 0,15, 0,TAU);
		Vec4 p = sph2cart(r*0.85f, TAU/4-0.25f, phi);
		vertex(vec4(p.z, p.y, p.x, 1.0f));
	}
	endLine();

	// laser - connections
	for(int i = 0; i < 15; i++) {
		float phi = map(i, 0,15, 0,TAU);
		v.pos = sph2cart(r, TAU/4-0.37f, phi);
		u.pos = sph2cart(r*0.85f, TAU/4-0.25f, phi);
		v.pos = vec4(v.pos.z, v.pos.y, v.pos.x, 1.0f);
		u.pos = vec4(u.pos.z, u.pos.y, u.pos.x, 1.0f);
		drawLine(v, u);
	}
}

void
setup(void)
{
	Serial.begin(9600);

	tft.initR(DISPTYPE);
	tft.setRotation(3);

	Serial.println(F("Initialized"));

	tft.fillScreen(BLACK);

	framebuffer = (u16*)malloc(WIDTH*HEIGHT*2);
	zbuffer = (u32*)malloc(WIDTH*HEIGHT*4);
	initDraw(framebuffer, zbuffer, WIDTH, HEIGHT);

	Model = m4ident();
	View = m4ident();
	Proj = m4persp(70.0f, ASPECT, 0.1f, 100.0f);
	updatemat();
}

void
scene_deathstar(void)
{
	float t = millis()/1000.0f;
	t += 4.0f;

	t *= -0.6f;
	Proj = m4ortho(2.0f, 2.0f, 0.1f, 10.0f);
	View = m4lookat(vec4(cos(t), sin(t), 0.2f, 1.0f),
	                vec4(0.0f, 0.0f, 0.0f, 1.0f),
	                vec4(0.0f, 0.0f, 1.0f, 0.0f));
	View = m4invOrtho(View);
	Model = m4ident();

	rs::srDepthTestEnable = 0;
	rs::srAlphaTestEnable = 0;
	rs::srTexEnable = 0;
	rs::srAlphaBlendEnable = 0;

	clear(1, 0);
	deathstar();
//	axes();
}

void
scene_cube(int fill)
{
	float t = millis()/1000.0f;

	Proj = m4persp(70.0f, ASPECT, 0.1f, 100.0f);
//	Proj = m4persp(70.0f, ASPECT, 2.0f, 10.0f);
	float d = 2.5f;
	View = m4lookat(vec4(d*cos(t), d*sin(t), 1.0f, 1.0f),
	                vec4(0.0f, 0.0f, 0.0f, 1.0f),
	                vec4(0.0f, 0.0f, 1.0f, 0.0f));
	View = m4invOrtho(View);
	Model = m4ident();

	rs::srDepthTestEnable = 1;
	rs::srAlphaTestEnable = 0;
	rs::srTexEnable = 0;
	rs::srAlphaBlendEnable = 0;

	t *= 0.6f;
//t = 0;
	Quat rx = qexp(quat(0,t,0,0));
	Quat ry = qexp(quat(0,0,t,0));
	Quat rz = qexp(quat(0,0,0,t));
	Model = quat2mat4(qmul(rx, qmul(ry, rz)));

	clear(1, 0);
	if(fill) {
		filled_cube();
	} else {
		cube();
		axes();
	}
}

void
loop(void)
{
	int sw = (millis()/10000) % 4;

//sw = 3;
	switch(sw) {
	case 0:
		scene_cube(0);
		break;
	case 1:
		scene_deathstar();
		break;
	case 2:
		scene_cube(0);
		break;
	case 3:
		scene_cube(1);
		break;
	}
	sendfb();
}
