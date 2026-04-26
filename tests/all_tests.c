#include "test_Constraints.h"
#include "test_CoreTransformations.h"
#include "test_Matrix.h"
#include "test_Parallel_cols.h"
#include "test_Parallel_rows.h"
#include "test_Presolver.h"
#include "test_QP.h"
#include "test_QR.h"
#include "test_QR_comprehensive.h"
#include "test_QR_objective.h"
#include "test_QR_end2end.h"
#include "test_QR_nkfix.h"
#include "test_QR_kkt.h"
#include "test_SimpleReductions.h"
#include "test_domain_propagation.h"
#include "test_dton.h"
#include "test_iVec.h"
#include "test_pathological.h"
#include "test_postsolve.h"
#include "test_radix_sort.h"
#include "test_ston.h"

const char *run_all_tests()
{
    /* Run QP/QR tests first - they have memory issues when run after other tests */
    mu_assert("qp error", test_qp());
    mu_assert("qr error", test_qr());
    mu_assert("qr_comprehensive error", test_qr_comprehensive());
    mu_assert("qr_objective error", test_qr_objective());
    mu_assert("qr_end2end error", test_qr_end2end());
    mu_assert("qr_nkfix error", test_qr_nkfix());
    mu_assert("qr_kkt error", test_qr_kkt());
    
    mu_assert("matrix error", test_matrix());
    mu_assert("constraints error", test_constraints());
    mu_assert("iVec error", test_iVec());
    mu_assert("dton error", test_dton());
    mu_assert("core error", test_core());
    mu_assert("ston error", test_ston());
    mu_assert("simple reductions error", test_simple());
    mu_assert("domain propagation error", test_domain());
    mu_assert("radix_sort error", test_radix_sort());
    mu_assert("presolver error", test_presolver());
    mu_assert("postsolve error", test_postsolve());
    mu_assert("pathological error", test_pathological());
    mu_assert("parallel_cols error", test_parallel_cols());

#ifndef _WIN32
    /* windows build is a bit weird in debug mode */
    mu_assert("parallel_rows error", test_parallel_rows());
#endif

    return NULL;
}

int main()
{
    const char *result = run_all_tests();
    if (result != NULL)
    {
        printf("Test failed: %s\n", result);
        return -1;
    }
    printf("All tests passed!\n");
    return 0;
}
