typedef struct { float x, y, z, w; } Vec4;
inline Vec4 vec4(float x, float y, float z, float w) { return (Vec4){x, y, z, w}; }
inline Vec4 v4neg(Vec4 v) { return vec4(-v.x, -v.y, -v.z, -v.w); }
inline Vec4 v4add(Vec4 u,Vec4 v) { return vec4(u.x+v.x, u.y+v.y, u.z+v.z, u.w+v.w); }
inline Vec4 v4sub(Vec4 u,Vec4 v) { return vec4(u.x-v.x, u.y-v.y, u.z-v.z, u.w-v.w); }
inline Vec4 v4scale(float s,Vec4 v) { return vec4(s*v.x, s*v.y, s*v.z, s*v.w); }
inline float v4dot(Vec4 a, Vec4 b) { return a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w; }
inline float v4normsq(Vec4 v) { return v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w; }
inline float v4norm(Vec4 v) { return sqrt(v4normsq(v)); }
inline Vec4 v4normalized(Vec4 v) { return v4scale(1/v4norm(v), v); }




typedef struct { Vec4 x, y, z, w; } Mat4;

inline Vec4 v4xform(Mat4 *m, Vec4 v) {
	return v4add(v4scale(v.x, m->x),
	       v4add(v4scale(v.y, m->y),
	       v4add(v4scale(v.z, m->z),
	             v4scale(v.w, m->w))));
}

inline Mat4 m4ident(void) {
	return (Mat4){
		{ 1.0f, 0.0f, 0.0f, 0.0f },
		{ 0.0f, 1.0f, 0.0f, 0.0f },
		{ 0.0f, 0.0f, 1.0f, 0.0f },
		{ 0.0f, 0.0f, 0.0f, 1.0f }
	};
}

inline Mat4 m4transposeP(Mat4 *m) {
	return (Mat4){
		{ m->x.x, m->y.x, m->z.x, m->w.x },
		{ m->x.y, m->y.y, m->z.y, m->w.y },
		{ m->x.z, m->y.z, m->z.z, m->w.z },
		{ m->x.w, m->y.w, m->z.w, m->w.w }
	};
}
inline Mat4 m4transpose(Mat4 m) { return m4transposeP(&m); }


inline Mat4 m4mulP(Mat4 *a, Mat4 *b) {
	return (Mat4){
		v4xform(a, b->x),
		v4xform(a, b->y),
		v4xform(a, b->z),
		v4xform(a, b->w)
	};
}
inline Mat4 m4mul(Mat4 a, Mat4 b) { return m4mulP(&a, &b); }

inline Mat4 m4invOrthoP(Mat4 *m) {
	return (Mat4){
		{ m->x.x, m->y.x, m->z.x, 0.0f },
		{ m->x.y, m->y.y, m->z.y, 0.0f },
		{ m->x.z, m->y.z, m->z.z, 0.0f },
		{ -v4dot(m->x, m->w),
		  -v4dot(m->y, m->w),
		  -v4dot(m->z, m->w),
		  1.0f }
	};
}
inline Mat4 m4invOrtho(Mat4 m) { return m4invOrthoP(&m); }


inline Mat4 m4persp(float fov, float aspect, float n, float f) {
	fov = fov/360.0f*6.28318f;
	float w = tan(fov/2.0f);
	float h = w/aspect;
	w = 1.0f/w;
	h = 1.0f/h;
#if 0
	// map [-n,-f]/w -> [0,1]
	float a = -f/(f-n);
	float b = -f*n/(f-n);
#else
	// map [-n,-f]/w -> [-1,1]
	float a = (n+f)/(n-f);
	float b = 2*n*f/(n-f);
#endif
	return (Mat4){
		{    w, 0.0f, 0.0f,  0.0f },
		{ 0.0f,    h, 0.0f,  0.0f },
		{ 0.0f, 0.0f,    a, -1.0f },
		{ 0.0f, 0.0f,    b,  0.0f }
	};
}

inline Mat4 m4ortho(float w, float h, float n, float f) {
	w = 2.0f/w;
	h = 2.0f/h;
	// map [-n,-f] -> [-1,1]
	float a = 2/(n-f);
	float b = (f+n)/(n-f);
	return (Mat4){
		{    w, 0.0f, 0.0f,  0.0f },
		{ 0.0f,    h, 0.0f,  0.0f },
		{ 0.0f, 0.0f,    a,  0.0f },
		{ 0.0f, 0.0f,    b,  1.0f }
	};
}

inline Mat4 m4scale(float x, float y, float z) {
	return (Mat4){
		{    x, 0.0f, 0.0f, 0.0f },
		{ 0.0f,    y, 0.0f, 0.0f },
		{ 0.0f, 0.0f,    z, 0.0f },
		{ 0.0f, 0.0f, 0.0f, 1.0f }
	};
}

// ugly hack :/
inline Vec4 v4cross3(Vec4 a, Vec4 b) {
	return vec4(
		a.y*b.z - a.z*b.y,
		a.z*b.x - a.x*b.z,
		a.x*b.y - a.y*b.x,
		0.0f);
}

inline Mat4 m4lookat(Vec4 pos, Vec4 target, Vec4 up) {
	Vec4 z = v4normalized(v4sub(pos, target));
	Vec4 x = v4normalized(v4cross3(up, z));
	Vec4 y = v4normalized(v4cross3(z, x));
	return (Mat4){ x, y, z, pos };
}
