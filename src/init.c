#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void Rborder(void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP createKDtree(SEXP, SEXP, SEXP);
extern SEXP RallRead(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Rarea(SEXP, SEXP);
extern SEXP Rballpivoting(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Rbarycenter(SEXP);
extern SEXP RBox(SEXP, SEXP);
extern SEXP Rclean(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RclosestKD(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Rclost(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RCone(SEXP, SEXP, SEXP, SEXP);
extern SEXP Rcurvature(SEXP, SEXP);
extern SEXP RDodecahedron(SEXP);
extern SEXP RgetEdge(SEXP, SEXP, SEXP);
extern SEXP RHexahedron(SEXP);
extern SEXP RIcosahedron(SEXP);
extern SEXP Rintersect(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Risolated(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Rkdtree(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Rkmeans(SEXP, SEXP, SEXP, SEXP);
extern SEXP RMarchC(SEXP, SEXP);
extern SEXP Rmeshres(SEXP, SEXP);
extern SEXP Rmeshvol(SEXP);
extern SEXP RmeshXPtr(SEXP);
extern SEXP Rmetro(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ROctahedron(SEXP);
extern SEXP ROneRing(SEXP, SEXP, SEXP);
extern SEXP RMeshWrite(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RQEdecim(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Rsample(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RsearchKDtree(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RsearchKDtreeForClosestPoints(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Rsmooth(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RSphere(SEXP, SEXP);
extern SEXP RSphericalCap(SEXP, SEXP, SEXP);
extern SEXP RSquare(SEXP);
// extern SEXP RSTLWrite(SEXP, SEXP, SEXP, SEXP);
extern SEXP Rsubdivision(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RTetrahedron(SEXP);
extern SEXP RuniformResampling(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RupdateNormals(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RVFadj(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"Rborder", (DL_FUNC) &Rborder, 6},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"createKDtree",                  (DL_FUNC) &createKDtree,                   3},
    {"RallRead",                      (DL_FUNC) &RallRead,                       6},
    {"Rarea",                         (DL_FUNC) &Rarea,                          2},
    {"Rballpivoting",                 (DL_FUNC) &Rballpivoting,                  5},
    {"Rbarycenter",                   (DL_FUNC) &Rbarycenter,                    1},
    {"RBox",                          (DL_FUNC) &RBox,                           2},
    {"Rclean",                        (DL_FUNC) &Rclean,                         5},
    {"RclosestKD",                    (DL_FUNC) &RclosestKD,                    13},
    {"Rclost",                        (DL_FUNC) &Rclost,                         8},
    {"RCone",                         (DL_FUNC) &RCone,                          4},
    {"Rcurvature",                    (DL_FUNC) &Rcurvature,                     2},
    {"RDodecahedron",                 (DL_FUNC) &RDodecahedron,                  1},
    {"RgetEdge",                      (DL_FUNC) &RgetEdge,                       3},
    {"RHexahedron",                   (DL_FUNC) &RHexahedron,                    1},
    {"RIcosahedron",                  (DL_FUNC) &RIcosahedron,                   1},
    {"Rintersect",                    (DL_FUNC) &Rintersect,                     8},
    {"Risolated",                     (DL_FUNC) &Risolated,                      6},
    {"Rkdtree",                       (DL_FUNC) &Rkdtree,                        6},
    {"Rkmeans",                       (DL_FUNC) &Rkmeans,                        4},
    {"RMarchC",                       (DL_FUNC) &RMarchC,                        2},
    {"Rmeshres",                      (DL_FUNC) &Rmeshres,                       2},
    {"Rmeshvol",                      (DL_FUNC) &Rmeshvol,                       1},
    {"RmeshXPtr",                     (DL_FUNC) &RmeshXPtr,                      1},
    {"Rmetro",                        (DL_FUNC) &Rmetro,                        14},
    {"ROctahedron",                   (DL_FUNC) &ROctahedron,                    1},
    {"ROneRing",                      (DL_FUNC) &ROneRing,                       3},
    {"RMeshWrite",                    (DL_FUNC) &RMeshWrite,                     8},
    {"RQEdecim",                      (DL_FUNC) &RQEdecim,                       5},
    {"Rsample",                       (DL_FUNC) &Rsample,                        5},
    {"RsearchKDtree",                 (DL_FUNC) &RsearchKDtree,                  5},
    {"RsearchKDtreeForClosestPoints", (DL_FUNC) &RsearchKDtreeForClosestPoints, 12},
    {"Rsmooth",                       (DL_FUNC) &Rsmooth,                        7},
    {"RSphere",                       (DL_FUNC) &RSphere,                        2},
    {"RSphericalCap",                 (DL_FUNC) &RSphericalCap,                  3},
    {"RSquare",                       (DL_FUNC) &RSquare,                        1},
    //    {"RSTLWrite",                     (DL_FUNC) &RSTLWrite,                      4},
    {"Rsubdivision",                  (DL_FUNC) &Rsubdivision,                   6},
    {"RTetrahedron",                  (DL_FUNC) &RTetrahedron,                   1},
    {"RuniformResampling",            (DL_FUNC) &RuniformResampling,             9},
    {"RupdateNormals",                (DL_FUNC) &RupdateNormals,                 5},
    {"RVFadj",                        (DL_FUNC) &RVFadj,                         2},
    {NULL, NULL, 0}
};

void R_init_Rvcg(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
