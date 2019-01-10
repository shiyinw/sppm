//#include "vec.hpp"
//#include <iostream>
//using namespace std;
//int M=100,N=100;
//Vec3f p[(M+1)*(N+1)];                          // 存储控制点
//std::vector<float> points(0);                  // 存储点
//std::vector<int> meshes(0);                    // 存储面
//int nu = 300, nv = 300;                         // 自定义采样密度
//
//for (int i = 0; i <= nu; i++) {
//    for (int j = 0; j <= nv; j++) {
//        float u = (float)i/nu, v = float(j)/nv;
//        points.push_back(Curve(u, v, p));
//        if (i != 0 && j != 0)
//            meshes.push_back(int4(.., .., .., ..)); // index从0开始
//    }
//}
