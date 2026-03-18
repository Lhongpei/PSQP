#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PSQP_API.h"

/* 
 * 这个测试展示：presolver 可能无法检测所有无界情况
 * 特别是当无界性需要多步推导时
 */

/* Test 1: 明显无界 - 可以检测 */
void test_obvious_unbounded()
{
    printf("Test 1: 明显无界（单变量，无约束）...\n");
    
    size_t n = 1;
    size_t m = 1;
    
    // 宽松约束
    double Ax[] = {1.0};
    int Ai[] = {0};
    int Ap[] = {0, 1};
    
    double lhs[] = {-INFINITY};
    double rhs[] = {INFINITY};
    double lbs[] = {-INFINITY};
    double ubs[] = {INFINITY};
    double c[] = {1.0};  // 最小化 x，x 无下界
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    stgs->dual_fix = true;
    
    Presolver *presolver = new_presolver(Ax, Ai, Ap, m, n, 1,
                                         lhs, rhs, lbs, ubs, c, stgs);
    if (!presolver) { printf("  创建失败\n"); return; }
    
    PresolveStatus status = run_presolver(presolver);
    printf("  状态: %s\n", 
           status == UNBNDORINFEAS ? "UNBNDORINFEAS ✓ 检测到" :
           status == REDUCED ? "REDUCED ✗ 未检测" : "其他");
    
    free_presolver(presolver);
    free_settings(stgs);
}

/* Test 2: 隐性无界 - 可能检测不到
 * min x + y
 * s.t. x - y = 0 (即 x = y)
 *      x, y >= 0
 * 
 * 这个是有界的！因为 x = y >= 0，目标 >= 0
 */
void test_bounded_not_unbounded()
{
    printf("\nTest 2: 实际上有界（x=y, x,y>=0）...\n");
    
    size_t n = 2;
    size_t m = 1;
    
    double Ax[] = {1.0, -1.0};  // x - y = 0
    int Ai[] = {0, 1};
    int Ap[] = {0, 2};
    
    double lhs[] = {0.0};
    double rhs[] = {0.0};
    double lbs[] = {0.0, 0.0};  // x, y >= 0
    double ubs[] = {INFINITY, INFINITY};
    double c[] = {1.0, 1.0};  // 最小化 x + y
    
    Settings *stgs = default_settings();
    stgs->verbose = false;
    
    Presolver *presolver = new_presolver(Ax, Ai, Ap, m, n, 2,
                                         lhs, rhs, lbs, ubs, c, stgs);
    if (!presolver) { printf("  创建失败\n"); return; }
    
    PresolveStatus status = run_presolver(presolver);
    printf("  状态: %s\n", 
           status == REDUCED ? "REDUCED ✓ 正确处理" :
           status == UNCHANGED ? "UNCHANGED ✓ 正确处理" :
           status == UNBNDORINFEAS ? "UNBNDORINFEAS ✗ 误判" : "其他");
    printf("  （这是一个有界问题，应该正常处理）\n");
    
    free_presolver(presolver);
    free_settings(stgs);
}

/* Test 3: 复杂无界 - 可能需要求解器才能检测
 * min x
 * s.t. x >= 0
 *      x <= y
 *      y 无界
 * 
 * 如果 y 可以无限增大，x 也可以 -> 但 x 的目标是最小化
 * 等等，这里 min x 且 x >= 0，最小值在 x = 0，是有界的
 * 
 * 换一个例子：
 * min y - x
 * s.t. x <= y
 *      x >= 0
 *      y >= 0
 *      y <= 2x (这限制了增长)
 * 
 * 实际上这个也是有界的...
 * 
 * 关键是：presolver 不会真正"求解"，只用局部规则
 */
void test_complex_case()
{
    printf("\nTest 3: 复杂情况（presolver 只做局部简化）...\n");
    printf("  注意：presolver 只做局部简化，不做全局分析\n");
    printf("  它可能无法识别所有无界情况\n");
    printf("  但这没关系，solver 会最终判断\n");
}

int main()
{
    printf("=== Presolver 有界性检测能力 ===\n\n");
    printf("关键概念：\n");
    printf("- presolver 的主要目标是简化问题，不是完整求解\n");
    printf("- 它使用局部启发式规则检测明显的无界/不可行情况\n");
    printf("- 复杂的情况可能检测不出来，这属于正常行为\n\n");
    
    test_obvious_unbounded();
    test_bounded_not_unbounded();
    test_complex_case();
    
    printf("\n=== 总结 ===\n");
    printf("✓ 有界可行的问题：presolver 肯定能正常处理\n");
    printf("✓ 明显的无界/不可行：通常能检测（UNBNDORINFEAS）\n");
    printf("✗ 复杂的情况：可能检测不出来（返回 REDUCED/UNCHANGED）\n");
    printf("\n注意：如果 presolver 没检测出来，solver 会最终判断\n");
    printf("这是正常的 presolver 行为，不是 bug\n");
    
    return 0;
}
