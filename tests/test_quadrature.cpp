//
// Created by Finn Eggers on 21.10.24.
//
#include "../src/math/quadrature.h"
#include <gtest/gtest.h>

using namespace fem::math;
using namespace fem::quadrature;

using Function       = std::function<Precision(Precision, Precision, Precision)>;
using FunctionList   = std::vector<Function>;
using FunctionResultList = std::vector<Precision>;
using DomainList     = std::vector<Domain>;

class QuadratureTest : public ::testing::Test {
    protected:
    std::map<Order, FunctionList> functions;
    std::map<Domain, std::map<Order, FunctionResultList>> domain_results;
    DomainList domain_list;

    void SetUp() override {
        functions = {
            {ORDER_CONSTANT,{
                [](double r, double s, double t) { return 1.0; },
                [](double r, double s, double t) { return 2.0; },
                [](double r, double s, double t) { return -1.0; },
                [](double r, double s, double t) { return 0.5; }}},

            {ORDER_LINEAR,{
                [](double r, double s, double t) { return r; },
                [](double r, double s, double t) { return s; },
                [](double r, double s, double t) { return t; },
                [](double r, double s, double t) { return r + s + t; }}},

            // Quadratic/bilinear functions
            {ORDER_QUADRATIC,{
                [](double r, double s, double t) { return r * r - r * s - t * s + 3.0 * t * t; },
                [](double r, double s, double t) { return r * s + r * t + s * t; },
                [](double r, double s, double t) { return r * r + s * s + t * t; },
                [](double r, double s, double t) { return 2.0 * r * r - 3.0 * s * s + 4.0 * t * s - r * s; }}},

            {ORDER_CUBIC,{
                [](double r, double s, double t) { return r * r * r + s * s * s + t * t * t + r * s; },
                [](double r, double s, double t) { return r * r * r + s * s * t + t * t * r + s; },
                [](double r, double s, double t) { return r * r * r - 2.0 * s * s * s + 3.0 * r * s * t + r; },
                [](double r, double s, double t) { return r * r * s - s * t * t + 4.0 * r * s * t + t;}}},

            {ORDER_QUARTIC,{
                [](double r, double s, double t) { return r * r * r * r + s * s * s * s + t * t * t * t; },
                [](double r, double s, double t) { return r * r * r * r - s * s * s * s + 2.0 * r * r * s * s + t * t * t * t; },
                [](double r, double s, double t) { return r * r * s * s + s * s * t * t + t * t * r * r; },
                [](double r, double s, double t) { return r * r * r * r - 2.0 * r * r * s * s + s * s * s * s - 3.0 * r * s * t * t; }}},

            {ORDER_QUINTIC,{
                [](double r, double s, double t) { return r * r * r * r * r + s * s * s * s * s + t * t * t * t * t + r * r + s * s; },
                [](double r, double s, double t) { return r * r * r * s * s + s * s * s * t * t + t * t * t * r * r + r * r + s * s; },
                [](double r, double s, double t) { return r * r * r * r * r - 3.0 * r * r * r * s * s + t * t * r * s * s + r * r + s * s; },
                [](double r, double s, double t) { return r * r * s * s * s + s * s * s * t * t - t * t * t * r * r + r + s; }}},
        };

        domain_results = {
            {DOMAIN_ISO_LINE_A, {
                                    {ORDER_CONSTANT, {2.00000000000000, 4.00000000000000, -2.00000000000000, 1.00000000000000}},
                                    {ORDER_LINEAR,   {0.0, 0.0, 0.0, 0.0}},
                                    {ORDER_QUADRATIC, {2.0 / 3.0, 0.0, 2.0 / 3.0, 1.33333333333333}},
                                    {ORDER_CUBIC,    {0.0, 0.0, 0.0, 0.0}},
                                    {ORDER_QUARTIC,  {2.0 / 5.0, 2.0 / 5.0, 0.0, 2.0 / 5.0}},
                                    {ORDER_QUINTIC,  {2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, 0.0}}
                                }},
            {DOMAIN_ISO_LINE_B, {
                                    {ORDER_CONSTANT, {1.00000000000000, 2.00000000000000, -1.00000000000000, 0.500000000000000}},
                                    {ORDER_LINEAR,   {1.0 / 2.0, 0.0, 0.0, 1.0 / 2.0}},
                                    {ORDER_QUADRATIC, {1.0 / 3.0, 0.0, 1.0 / 3.0, 0.666666666666667}},
                                    {ORDER_CUBIC,    {1.0 / 4.0, 1.0 / 4.0, 3.0 / 4.0, 0.0}},
                                    {ORDER_QUARTIC,  {1.0 / 5.0, 1.0 / 5.0, 0.0, 1.0 / 5.0}},
                                    {ORDER_QUINTIC,  {1.0 / 2.0, 1.0 / 3.0, 1.0 / 2.0, 1.0 / 2.0}}
                                }},
            {DOMAIN_ISO_TRI, {
                                 {ORDER_CONSTANT, {0.500000000000000, 1.00000000000000, -0.500000000000000, 0.250000000000000}},
                                 {ORDER_LINEAR,   {1.0 / 6.0, 1.0 / 6.0, 0.0, 1.0 / 3.0}},
                                 {ORDER_QUADRATIC, {1.0 / 24.0, 1.0 / 24.0, 1.0 / 6.0, -0.125000000000000}},
                                 {ORDER_CUBIC,    {17.0 / 120.0, 13.0 / 60.0, 0.116666666666667, 1.0 / 60.0}},
                                 {ORDER_QUARTIC,  {1.0 / 15.0, 0.0111111111111111, 1.0 / 180.0, 0.0555555555555556}},
                                 {ORDER_QUINTIC,  {3.0 / 14.0, 71.0 / 420.0, 0.183333333333334, 47.0 / 140.0}}
                             }},
            {DOMAIN_ISO_QUAD, {
                                  {ORDER_CONSTANT, {4.00000000000000, 8.00000000000000, -4.00000000000000, 2.00000000000000}},
                                  {ORDER_LINEAR,   {0.0, 0.0, 0.0, 0.0}},
                                  {ORDER_QUADRATIC, {4.0 / 3.0, 0.0, 8.0 / 3.0, -1.33333333333333}},
                                  {ORDER_CUBIC,    {0.0, 0.0, 0.0, 0.0}},
                                  {ORDER_QUARTIC,  {8.0 / 5.0, 0.888888888888889, 4.0 / 9.0, 0.711111111111111}},
                                  {ORDER_QUINTIC,  {8.0 / 3.0, 8.0 / 3.0, 2.66666666666667, 0.0}}
                              }},
            {DOMAIN_ISO_HEX, {
                                 {ORDER_CONSTANT, {8.00000000000000, 16.0000000000000, -8.00000000000000, 4.00000000000000}},
                                 {ORDER_LINEAR,   {0.0, 0.0, 0.0, 0.0}},
                                 {ORDER_QUADRATIC, {10.6666666666667, 0.0, 8.0, -2.66666666666667}},
                                 {ORDER_CUBIC,    {0.0, 0.0, 0.0, 0.0}},
                                 {ORDER_QUARTIC,  {24.0 / 5.0, 3.37777777777778, 8.0 / 3.0, 1.42222222222222}},
                                 {ORDER_QUINTIC,  {16.0 / 3.0, 16.0 / 3.0, 5.33333333333333, 0.0}}
                             }},
            {DOMAIN_ISO_TET, {
                                 {ORDER_CONSTANT, {0.166666666666667, 0.333333333333333, -0.166666666666667, 0.0833333333333333}},
                                 {ORDER_LINEAR,   {1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 8.0}},
                                 {ORDER_QUADRATIC, {0.0500000000000003, 1.0 / 40.0, 1.0 / 20.0, 0.00833333333333308}},
                                 {ORDER_CUBIC,    {1.0 / 30.0, 1.0 / 18.0, 0.0375000000000001, 0.0472222222222219}},
                                 {ORDER_QUARTIC,  {1.0 / 70.0, 0.00634920634920638, 1.0 / 420.0, 0.00674603174603161}},
                                 {ORDER_QUINTIC,  {71.0 / 1680.0, 23.0 / 672.0, 0.0355158730158730, 281.0 / 3360.0}}
                             }},
            {DOMAIN_ISO_WEDGE, {
                                   {ORDER_CONSTANT, {1.00000000000000, 2.00000000000000, -1.00000000000000, 0.500000000000000}},
                                   {ORDER_LINEAR,   {1.0 / 3.0, 1.0 / 3.0, 0.0, 2.0 / 3.0}},
                                   {ORDER_QUADRATIC, {1.08333333333333, 1.0 / 12.0, 2.0 / 3.0, -0.250000000000000}},
                                   {ORDER_CUBIC,    {17.0 / 60.0, 49.0 / 90.0, 0.233333333333333, -0.0777777777777778}},
                                   {ORDER_QUARTIC,  {1.0 / 3.0, 0.222222222222222, 11.0 / 90.0, 0.0277777777777779}},
                                   {ORDER_QUINTIC,  {3.0 / 7.0, 13.0 / 35.0, 0.377777777777778, 74.0 / 105.0}}
                               }},
            {DOMAIN_ISO_PYRAMID, {
                                     {ORDER_CONSTANT, {1.33333333333333, 2.66666666666667, -1.33333333333333, 0.666666666666667}},
                                     {ORDER_LINEAR,   {0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0}},
                                     {ORDER_QUADRATIC, {0.666666666666669, 0.0, 2.0 / 3.0, -0.266666666666666}},
                                     {ORDER_CUBIC,    {1.0 / 15.0, 2.0 / 45.0, 0.0, 0.333333333333333}},
                                     {ORDER_QUARTIC,  {4.0 / 15.0, 0.165079365079364, 4.0 / 45.0, 0.101587301587303}},
                                     {ORDER_QUINTIC,  {39.0 / 70.0, 113.0 / 210.0, 0.533333333333334, -1.0 / 210.0}}
                                 }}
        };


        domain_list = {
            DOMAIN_ISO_LINE_A,
            DOMAIN_ISO_LINE_B,
            DOMAIN_ISO_TRI,
            DOMAIN_ISO_QUAD,
            DOMAIN_ISO_HEX,
            DOMAIN_ISO_TET,
            DOMAIN_ISO_WEDGE,
            DOMAIN_ISO_PYRAMID
        };
    }
};
#define TEST_DOMAIN_ORDER_FUNCTION(domain, order, func_order, idx) \
    TEST_F(QuadratureTest, _##domain##_##order##_FUNCTION_##func_order##_##idx) { \
        try {                                                     \
            auto true_result = domain_results[domain][func_order][idx]; \
            auto test_result = Quadrature(domain, order).integrate(functions[func_order][idx]);\
            EXPECT_NEAR(test_result, true_result, 1e-10);         \
        } catch (std::runtime_error& e) {                         \
        }                                                         \
    }


#define TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, order, func_order)\
    TEST_DOMAIN_ORDER_FUNCTION(domain, order, func_order, 0); \
    TEST_DOMAIN_ORDER_FUNCTION(domain, order, func_order, 1);\
    TEST_DOMAIN_ORDER_FUNCTION(domain, order, func_order, 2);\
    TEST_DOMAIN_ORDER_FUNCTION(domain, order, func_order, 3);\

#define TEST_DOMAIN_ORDER_CONSTANT(domain)\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_CONSTANT, ORDER_CONSTANT);

#define TEST_DOMAIN_ORDER_LINEAR(domain)\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_LINEAR, ORDER_CONSTANT);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_LINEAR, ORDER_LINEAR);\

#define TEST_DOMAIN_ORDER_QUADRATIC(domain)\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_QUADRATIC, ORDER_CONSTANT);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_QUADRATIC, ORDER_LINEAR);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_QUADRATIC, ORDER_QUADRATIC);\

#define TEST_DOMAIN_ORDER_CUBIC(domain)\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_CUBIC, ORDER_CONSTANT);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_CUBIC, ORDER_LINEAR);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_CUBIC, ORDER_QUADRATIC);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_CUBIC, ORDER_CUBIC);\

#define TEST_DOMAIN_ORDER_QUARTIC(domain)\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_QUARTIC, ORDER_CONSTANT);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_QUARTIC, ORDER_LINEAR);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_QUARTIC, ORDER_QUADRATIC);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_QUARTIC, ORDER_CUBIC);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_QUARTIC, ORDER_QUARTIC);\

#define TEST_DOMAIN_ORDER_QUINTIC(domain)\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_QUINTIC, ORDER_CONSTANT);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_QUINTIC, ORDER_LINEAR);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_QUINTIC, ORDER_QUADRATIC);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_QUINTIC, ORDER_CUBIC);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_QUINTIC, ORDER_QUARTIC);\
    TEST_DOMAIN_ORDER_CONSTANT_FUNCTIONS(domain, ORDER_QUINTIC, ORDER_QUINTIC);

#define TEST_DOMAIN(domain)\
    TEST_DOMAIN_ORDER_CONSTANT(domain);\
    TEST_DOMAIN_ORDER_LINEAR(domain);\
    TEST_DOMAIN_ORDER_QUADRATIC(domain);\
    TEST_DOMAIN_ORDER_CUBIC(domain);\
    TEST_DOMAIN_ORDER_QUARTIC(domain);\
    TEST_DOMAIN_ORDER_QUINTIC(domain);

TEST_DOMAIN(DOMAIN_ISO_LINE_A)
TEST_DOMAIN(DOMAIN_ISO_LINE_B)
TEST_DOMAIN(DOMAIN_ISO_TRI)
TEST_DOMAIN(DOMAIN_ISO_QUAD)
TEST_DOMAIN(DOMAIN_ISO_HEX)
TEST_DOMAIN(DOMAIN_ISO_TET)
TEST_DOMAIN(DOMAIN_ISO_WEDGE)
TEST_DOMAIN(DOMAIN_ISO_PYRAMID)