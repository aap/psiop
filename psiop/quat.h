typedef struct { float w, x, y, z; } Quat;
inline Quat quat(float r, float x, float y, float z) { return (Quat){r, x, y, z}; }
inline Quat rquat(float r) { return (Quat){r, 0.0f, 0.0f, 0.0f}; }
inline Quat qvec(Quat q) { return quat(0, q.x, q.y, q.z); }
inline Quat qneg(Quat q) { return quat(-q.w, -q.x, -q.y, -q.z); }
inline Quat qadd(Quat p,Quat q) { return quat(p.w+q.w, p.x+q.x, p.y+q.y, p.z+q.z); }
inline Quat qsub(Quat p,Quat q) { return quat(p.w-q.w, p.x-q.x, p.y-q.y, p.z-q.z); }
inline Quat qscale(float s,Quat q) { return quat(s*q.w, s*q.x, s*q.y, s*q.z); }
inline Quat qconj(Quat q) { return quat(q.w, -q.x, -q.y, -q.z); }
inline Quat qmul(Quat p,Quat q) { return quat(
  p.w*q.w - p.x*q.x - p.y*q.y - p.z*q.z,
  p.w*q.x + p.x*q.w + p.y*q.z - p.z*q.y,
  p.w*q.y + p.y*q.w + p.z*q.x - p.x*q.z,
  p.w*q.z + p.z*q.w + p.x*q.y - p.y*q.x
); }
inline float qnormsq(Quat q) { return q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z; }
inline float qnorm(Quat q) { return sqrt(qnormsq(q)); }
inline Quat qinv(Quat q) { return qscale(1/qnormsq(q), qconj(q)); }
inline Quat qnormalized(Quat q) { return qscale(1/qnorm(q), q); }
inline Quat qsqrtnorm(Quat q) { return qnormalized(qadd(q, rquat(1.0f))); }	// q normalized before and after
inline Quat qfromto(Quat u, Quat v) { return qsqrtnorm(qmul(quat(0,v.x,v.y,v.z), quat(0,-u.x,-u.y,-u.z))); }
inline Quat qsandwich(Quat q, Quat v) { return qmul(q, qmul(v, qconj(q))); }
inline Quat qrsandwich(Quat q, Quat v) { return qmul(qconj(q), qmul(v, q)); }
inline Quat qrotate(Quat q, Quat v) { return qmul(q, qmul(qvec(v), qinv(q))); }
inline Quat qrrotate(Quat q, Quat v) { return qmul(qinv(q), qmul(qvec(v), q)); }
inline Quat qexp(Quat q) {
	Quat qv = qvec(q);
	float l = qnorm(qv);
	if(l == 0.0f)
		return rquat(exp(q.w));
	return qscale(exp(q.w), qadd(rquat(cos(l)), qscale(sin(l)/l, qv)));
}
inline Quat qlog(Quat q) {
	float l = qnorm(q);
	q = qscale(1.0f/l, q);
	float c = q.w; q.w = 0;
	float s = qnorm(q);
	if(s == 0.0f)
		return rquat(log(l));
	return qadd(rquat(log(l)), qscale(atan2(s,c)/s, q));
}
inline Quat qpow(Quat q, float a) { return qexp(qscale(a, qlog(q))); }
inline Quat qcross(Quat p,Quat q) { return qscale(0.5f, qsub(qmul(p,q), qmul(q,p))); }
