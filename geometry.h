template<typename Type>
class Vec3 {

public:
    Type x,y,z;

    Vec3() {}
    Vec3(Type X, Type Y, Type Z) : x(X), y(Y), z(Z) {}

    Vec3 operator + (Vec3 &V)
    {
        return Vec3(x + V.x, y + V.y, z + V.z); 
    }

    Type operator [] (int i)
    {
        return (&x)[i];
    }
};

template<typename Type>
class Vec2 {

};

template<typename Type>
class Mat44{

public:
    void inverse();

    void multiple(Vec3<Type> &Vec, Vec3<Type> &Result)
    {
    /* Vec * Mat44*/

    Type a, b, c, w;

    a = Vec[0] * m[0][0] + Vec[1] * m[1][0] + Vec[2] * m[2][0] + m[3][0];
    b = Vec[0] * m[0][1] + Vec[1] * m[1][1] + Vec[2] * m[2][1] + m[3][1];
    c = Vec[0] * m[0][2] + Vec[1] * m[1][2] + Vec[2] * m[2][2] + m[3][2];
    w = Vec[0] * m[0][3] + Vec[1] * m[1][3] + Vec[2] * m[2][3] + m[3][3];

    Result.x = a/w;
    Result.y = a/w;
    Result.z = a/w;
    }


    Type m[4][4];

};

typedef Vec3<float>  Vec3f;
typedef Vec2<int>    Vec2i;
typedef Mat44<float> Mat44f;

