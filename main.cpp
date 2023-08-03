#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
#include <set>
#include <math.h>

struct __attribute__((__packed__)) STL_tri {
    float norm[3];
    float p1[3];
    float p2[3];
    float p3[3];
    short int v;
};

class STLtoTAB {
public:
    float dx,dy;
    int nx, ny;
    std::vector< double > tab;
    STLtoTAB(): dx(0),dy(0),nx(0),ny(0) {

    }
    void add (const float& x,const float& y,const float& z) {
        int ix = round(x/dx);
        int iy = round(y/dy);
        if (fabs(x/dx - ix) > 1e-3 || fabs(y/dy - iy) > 1e-3) {
            printf("doesn't align: %g/%g, %g/%g\n",x,dx,y,dy);
            exit(-1);
        }
        if (ix >= nx || iy >= ny) {
            printf("out-of-bounds: %d/%d, %d/%d\n",ix,nx,iy,ny);
            printf("point: %g/%g, %g/%g\n",x,dx,y,dy);
            exit(-1);
        }
        double val = tab[ix+nx*iy];
        if (val != NAN) {
            if (fabs(val - z)/(val+z+1e-6) > 1e-3) {
                printf("conflict value at: %d, %d\n",ix,iy);
                printf("new: %g, old: %g\n",z, val);
                exit(-1);
            }
        }
        tab[ix+nx*iy] = z;
    }
    void load(const std::vector<STL_tri>& tri) {
        dx = tri[0].p2[0];
        dy = tri[0].p2[1];
        printf(" DX: %g DY: %g\n", dx,dy);
        float mx = 0, my = 0;
        for (auto t : tri) {
            if (t.p1[0] > mx) mx = t.p1[0];
            if (t.p2[0] > mx) mx = t.p2[0];
            if (t.p3[0] > mx) mx = t.p3[0];
            if (t.p1[1] > my) my = t.p1[1];
            if (t.p2[1] > my) my = t.p2[1];
            if (t.p3[1] > my) my = t.p3[1];
        }
        nx = mx / dx + 1;
        ny = my / dy + 1;
        printf(" NX: %d NY: %d\n", nx,ny);
        tab.resize(nx*ny, NAN);
        for (auto t : tri) {
            add(t.p1[0],t.p1[1],t.p1[2]);
            add(t.p2[0],t.p2[1],t.p2[2]);
            add(t.p3[0],t.p3[1],t.p3[2]);
        }
    }
};

class STL {
public:
    char header[80];
    std::vector<STL_tri> tri;
    STL () { }
    void load(std::string filename) {
        const char * fn = filename.c_str();  
        int ret,ntri;
        FILE *f = fopen(fn, "rb");
        if (f == NULL) {
            printf("'STL' element: %s doesn't exists or cannot be opened\n", fn);
            exit(-1);
        }
        ret = fread(header, 80, sizeof(char), f);
        ret = fread(&ntri, 1, sizeof(int), f);
        if (!strncmp(header, "solid", 5)){      // Checking if STL is binary. STL in ASCII begins with "solid"
            printf("'STL' element %s is not in binary format!\n", fn);
            exit(-1);
        }
        printf("Alloc...\n");
        tri.resize(ntri);
        printf("Reading...\n");
        ret = fread((void*) tri.data(), ntri, sizeof(STL_tri), f);
        printf("Closing...\n");
        fclose(f);
    }
};

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix read_stl_to_tab(std::string filename) {
    STL stl;
    stl.load(filename);
    STLtoTAB tab;
    tab.load(stl.tri);
    NumericMatrix ret( tab.nx , tab.ny , tab.tab.begin() );
    ret.attr("pixel") = NumericVector({tab.dx, tab.dy});
    return ret;
}
