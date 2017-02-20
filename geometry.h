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

    friend std::ostream& operator << (std::ostream &S, const Vec3<Type> &V)
    {
        return S << '[' << V.x << ' ' << V.y << ' ' << V.z << ']';
    }
};

template<typename Type>
class Vec2 {

public:
    Type x,y,z;

    friend std::ostream& operator << (std::ostream &S, const Vec2<Type> &V)
    {
        return S << '[' << V.x << ' ' << V.y << ']';
    }
};

template<typename Type>
class Mat44{

public:
    Mat44()
    {

    }

    Mat44(Type a, Type b, Type c, Type d, Type e, Type f, Type g, Type h, Type i, Type j, Type k, Type l, Type m, Type n, Type o, Type p)
    {
        x[0][0] = a;
        x[0][1] = b;
        x[0][2] = c;
        x[0][3] = d;
        x[1][0] = e;
        x[1][1] = f;
        x[1][2] = g;
        x[1][3] = h;
        x[2][0] = i;
        x[2][1] = j;
        x[2][2] = k;
        x[2][3] = l;
        x[3][0] = m;
        x[3][1] = n;
        x[3][2] = o;
        x[3][3] = p;
    }

    Mat44 inverse()
    {
        int i, j, k;
        Mat44 s;
        Mat44 t (*this);
        
        // Forward elimination
        for (i = 0; i < 3 ; i++) {
            int pivot = i;
            
            Type pivotsize = t[i][i];
            
            if (pivotsize < 0)
                pivotsize = -pivotsize;
                
                for (j = i + 1; j < 4; j++) {
                    Type tmp = t[j][i];
                    
                    if (tmp < 0)
                        tmp = -tmp;
                        
                        if (tmp > pivotsize) {
                            pivot = j;
                            pivotsize = tmp;
                        }
                }
            
            if (pivotsize == 0) {
                // Cannot invert singular matrix
                return Mat44();
            }
            
            if (pivot != i) {
                for (j = 0; j < 4; j++) {
                    Type tmp;
                    
                    tmp = t[i][j];
                    t[i][j] = t[pivot][j];
                    t[pivot][j] = tmp;
                    
                    tmp = s[i][j];
                    s[i][j] = s[pivot][j];
                    s[pivot][j] = tmp;
                }
            }
            
            for (j = i + 1; j < 4; j++) {
                Type f = t[j][i] / t[i][i];
                
                for (k = 0; k < 4; k++) {
                    t[j][k] -= f * t[i][k];
                    s[j][k] -= f * s[i][k];
                }
            }
        }
        
        // Backward substitution
        for (i = 3; i >= 0; --i) {
            Type f;
            
            if ((f = t[i][i]) == 0) {
                // Cannot invert singular matrix
                return Mat44();
            }
            
            for (j = 0; j < 4; j++) {
                t[i][j] /= f;
                s[i][j] /= f;
            }
            
            for (j = 0; j < i; j++) {
                f = t[j][i];
                
                for (k = 0; k < 4; k++) {
                    t[j][k] -= f * t[i][k];
                    s[j][k] -= f * s[i][k];
                }
            }
        }
        
        return s; 
    }

    void multiple(Vec3<Type> &Vec, Vec3<Type> &Result)
    {
        /* Vec * Mat44*/

        Type a, b, c, w;

        a = Vec[0] * x[0][0] + Vec[1] * x[1][0] + Vec[2] * x[2][0] + x[3][0];
        b = Vec[0] * x[0][1] + Vec[1] * x[1][1] + Vec[2] * x[2][1] + x[3][1];
        c = Vec[0] * x[0][2] + Vec[1] * x[1][2] + Vec[2] * x[2][2] + x[3][2];
        w = Vec[0] * x[0][3] + Vec[1] * x[1][3] + Vec[2] * x[2][3] + x[3][3];

        Result.x = a/w;
        Result.y = b/w;
        Result.z = c/w;
    }

    Type* operator[] (int i)
    {
        return x[i];
    }

    const Type* operator[] (int i) const
    {
        return x[i];
    }

    friend std::ostream& operator << (std::ostream &s, const Mat44<Type> &m)
    {
        std::ios_base::fmtflags oldFlags = s.flags();
        int width = 12; // total with of the displayed number
        s.precision(5); // control the number of displayed decimals
        s.setf (std::ios_base::fixed);
        
        s << "[" << std::setw (width) << m[0][0] <<
             " " << std::setw (width) << m[0][1] <<
             " " << std::setw (width) << m[0][2] <<
             " " << std::setw (width) << m[0][3] << "\n" <<
            
             " " << std::setw (width) << m[1][0] <<
             " " << std::setw (width) << m[1][1] <<
             " " << std::setw (width) << m[1][2] <<
             " " << std::setw (width) << m[1][3] << "\n" <<
            
             " " << std::setw (width) << m[2][0] <<
             " " << std::setw (width) << m[2][1] <<
             " " << std::setw (width) << m[2][2] <<
             " " << std::setw (width) << m[2][3] << "\n" <<
            
             " " << std::setw (width) << m[3][0] <<
             " " << std::setw (width) << m[3][1] <<
             " " << std::setw (width) << m[3][2] <<
             " " << std::setw (width) << m[3][3] << "]";
        
        s.flags (oldFlags);
        return s;
    }

    Type x[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

};

typedef Vec3<float>  Vec3f;
typedef Vec2<int>    Vec2i;
typedef Vec2<float>  Vec2f;
typedef Mat44<float> Mat44f;
