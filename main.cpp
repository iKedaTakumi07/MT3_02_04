#include <Novice.h>
#define _USE_MATH_DEFINES
#include <assert.h>
#include <cmath>
#include <cstdint>
#include <imgui.h>
#include <math.h>

const char kWindowTitle[] = "LE2C_03_イケダ_タクミ";

struct Vector3 {
    float x;
    float y;
    float z;
};
struct Matrix4x4 {
    float m[4][4];
};
struct Sphere {
    Vector3 center;
    float radius;
    int color;
};
struct Plane {
    Vector3 normal;
    float distance;
};
struct Segment {
    Vector3 origin;
    Vector3 diff;
};
struct Triangle {
    Vector3 vertices[3];
};

Vector3 add(const Vector3& v1, const Vector3& v2)
{
    Vector3 num;
    num.x = v1.x + v2.x;
    num.y = v1.y + v2.y;
    num.z = v1.z + v2.z;
    return num;
}
Vector3 Subtract(const Vector3& v1, const Vector3& v2)
{
    Vector3 num;
    num.x = v1.x - v2.x;
    num.y = v1.y - v2.y;
    num.z = v1.z - v2.z;
    return num;
}
float Dot(const Vector3& v1, const Vector3& v2)
{
    float num;
    return num = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
float Length(const Vector3& v)
{
    float num;
    return num = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2)
{
    Matrix4x4 num;
    num.m[0][0] = m1.m[0][0] * m2.m[0][0] + m1.m[0][1] * m2.m[1][0] + m1.m[0][2] * m2.m[2][0] + m1.m[0][3] * m2.m[3][0];
    num.m[0][1] = m1.m[0][0] * m2.m[0][1] + m1.m[0][1] * m2.m[1][1] + m1.m[0][2] * m2.m[2][1] + m1.m[0][3] * m2.m[3][1];
    num.m[0][2] = m1.m[0][0] * m2.m[0][2] + m1.m[0][1] * m2.m[1][2] + m1.m[0][2] * m2.m[2][2] + m1.m[0][3] * m2.m[3][2];
    num.m[0][3] = m1.m[0][0] * m2.m[0][3] + m1.m[0][1] * m2.m[1][3] + m1.m[0][2] * m2.m[2][3] + m1.m[0][3] * m2.m[3][3];

    num.m[1][0] = m1.m[1][0] * m2.m[0][0] + m1.m[1][1] * m2.m[1][0] + m1.m[1][2] * m2.m[2][0] + m1.m[1][3] * m2.m[3][0];
    num.m[1][1] = m1.m[1][0] * m2.m[0][1] + m1.m[1][1] * m2.m[1][1] + m1.m[1][2] * m2.m[2][1] + m1.m[1][3] * m2.m[3][1];
    num.m[1][2] = m1.m[1][0] * m2.m[0][2] + m1.m[1][1] * m2.m[1][2] + m1.m[1][2] * m2.m[2][2] + m1.m[1][3] * m2.m[3][2];
    num.m[1][3] = m1.m[1][0] * m2.m[0][3] + m1.m[1][1] * m2.m[1][3] + m1.m[1][2] * m2.m[2][3] + m1.m[1][3] * m2.m[3][3];

    num.m[2][0] = m1.m[2][0] * m2.m[0][0] + m1.m[2][1] * m2.m[1][0] + m1.m[2][2] * m2.m[2][0] + m1.m[2][3] * m2.m[3][0];
    num.m[2][1] = m1.m[2][0] * m2.m[0][1] + m1.m[2][1] * m2.m[1][1] + m1.m[2][2] * m2.m[2][1] + m1.m[2][3] * m2.m[3][1];
    num.m[2][2] = m1.m[2][0] * m2.m[0][2] + m1.m[2][1] * m2.m[1][2] + m1.m[2][2] * m2.m[2][2] + m1.m[2][3] * m2.m[3][2];
    num.m[2][3] = m1.m[2][0] * m2.m[0][3] + m1.m[2][1] * m2.m[1][3] + m1.m[2][2] * m2.m[2][3] + m1.m[2][3] * m2.m[3][3];

    num.m[3][0] = m1.m[3][0] * m2.m[0][0] + m1.m[3][1] * m2.m[1][0] + m1.m[3][2] * m2.m[2][0] + m1.m[3][3] * m2.m[3][0];
    num.m[3][1] = m1.m[3][0] * m2.m[0][1] + m1.m[3][1] * m2.m[1][1] + m1.m[3][2] * m2.m[2][1] + m1.m[3][3] * m2.m[3][1];
    num.m[3][2] = m1.m[3][0] * m2.m[0][2] + m1.m[3][1] * m2.m[1][2] + m1.m[3][2] * m2.m[2][2] + m1.m[3][3] * m2.m[3][2];
    num.m[3][3] = m1.m[3][0] * m2.m[0][3] + m1.m[3][1] * m2.m[1][3] + m1.m[3][2] * m2.m[2][3] + m1.m[3][3] * m2.m[3][3];

    return num;
}
Vector3 Multiply(const float& m2, const Vector3& m1)
{
    Vector3 num;
    num.x = m1.x * m2;
    num.y = m1.y * m2;
    num.z = m1.z * m2;

    return num;
}

Matrix4x4 Inverse(const Matrix4x4& m)
{
    float determinant;
    Matrix4x4 num;

    determinant = m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1] + m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2]
        - m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1] - m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3] - m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2]
        - m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1] - m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2]
        + m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1] + m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2]
        + m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1] + m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2]
        - m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1] - m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3] - m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2]
        - m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0] - m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0] - m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0]
        + m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0] + m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0] + m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0];

    if (determinant == 0.0f) {
        return m;
    };

    num.m[0][0] = (m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[1][2] * m.m[2][3] * m.m[3][1] + m.m[1][3] * m.m[2][1] * m.m[3][2] - m.m[1][3] * m.m[2][2] * m.m[3][1] - m.m[1][2] * m.m[2][1] * m.m[3][3] - m.m[1][1] * m.m[2][3] * m.m[3][2]) / determinant;
    num.m[0][1] = (-m.m[0][1] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[2][3] * m.m[3][1] - m.m[0][3] * m.m[2][1] * m.m[3][2] + m.m[0][3] * m.m[2][2] * m.m[3][1] + m.m[0][2] * m.m[2][1] * m.m[3][3] + m.m[0][1] * m.m[2][3] * m.m[3][2]) / determinant;
    num.m[0][2] = (m.m[0][1] * m.m[1][2] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[3][1] + m.m[0][3] * m.m[1][1] * m.m[3][2] - m.m[0][3] * m.m[1][2] * m.m[3][1] - m.m[0][2] * m.m[1][1] * m.m[3][3] - m.m[0][1] * m.m[1][3] * m.m[3][2]) / determinant;
    num.m[0][3] = (-m.m[0][1] * m.m[1][2] * m.m[2][3] - m.m[0][2] * m.m[1][3] * m.m[2][1] - m.m[0][3] * m.m[1][1] * m.m[2][2] + m.m[0][3] * m.m[1][2] * m.m[2][1] + m.m[0][2] * m.m[1][1] * m.m[2][3] + m.m[0][1] * m.m[1][3] * m.m[2][2]) / determinant;

    num.m[1][0] = (-m.m[1][0] * m.m[2][2] * m.m[3][3] - m.m[1][2] * m.m[2][3] * m.m[3][0] - m.m[1][3] * m.m[2][0] * m.m[3][2] + m.m[1][3] * m.m[2][2] * m.m[3][0] + m.m[1][2] * m.m[2][0] * m.m[3][3] + m.m[1][0] * m.m[2][3] * m.m[3][2]) / determinant;
    num.m[1][1] = (m.m[0][0] * m.m[2][2] * m.m[3][3] + m.m[0][2] * m.m[2][3] * m.m[3][0] + m.m[0][3] * m.m[2][0] * m.m[3][2] - m.m[0][3] * m.m[2][2] * m.m[3][0] - m.m[0][2] * m.m[2][0] * m.m[3][3] - m.m[0][0] * m.m[2][3] * m.m[3][2]) / determinant;
    num.m[1][2] = (-m.m[0][0] * m.m[1][2] * m.m[3][3] - m.m[0][2] * m.m[1][3] * m.m[3][0] - m.m[0][3] * m.m[1][0] * m.m[3][2] + m.m[0][3] * m.m[1][2] * m.m[3][0] + m.m[0][2] * m.m[1][0] * m.m[3][3] + m.m[0][0] * m.m[1][3] * m.m[3][2]) / determinant;
    num.m[1][3] = (m.m[0][0] * m.m[1][2] * m.m[2][3] + m.m[0][2] * m.m[1][3] * m.m[2][0] + m.m[0][3] * m.m[1][0] * m.m[2][2] - m.m[0][3] * m.m[1][2] * m.m[2][0] - m.m[0][2] * m.m[1][0] * m.m[2][3] - m.m[0][0] * m.m[1][3] * m.m[2][2]) / determinant;

    num.m[2][0] = (m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[1][1] * m.m[2][3] * m.m[3][0] + m.m[1][3] * m.m[2][0] * m.m[3][1] - m.m[1][3] * m.m[2][1] * m.m[3][0] - m.m[1][1] * m.m[2][0] * m.m[3][3] - m.m[1][0] * m.m[2][3] * m.m[3][1]) / determinant;
    num.m[2][1] = (-m.m[0][0] * m.m[2][1] * m.m[3][3] - m.m[0][1] * m.m[2][3] * m.m[3][0] - m.m[0][3] * m.m[2][0] * m.m[3][1] + m.m[0][3] * m.m[2][1] * m.m[3][0] + m.m[0][1] * m.m[2][0] * m.m[3][3] + m.m[0][0] * m.m[2][3] * m.m[3][1]) / determinant;
    num.m[2][2] = (m.m[0][0] * m.m[1][1] * m.m[3][3] + m.m[0][1] * m.m[1][3] * m.m[3][0] + m.m[0][3] * m.m[1][0] * m.m[3][1] - m.m[0][3] * m.m[1][1] * m.m[3][0] - m.m[0][1] * m.m[1][0] * m.m[3][3] - m.m[0][0] * m.m[1][3] * m.m[3][1]) / determinant;
    num.m[2][3] = (-m.m[0][0] * m.m[1][1] * m.m[2][3] - m.m[0][1] * m.m[1][3] * m.m[2][0] - m.m[0][3] * m.m[1][0] * m.m[2][1] + m.m[0][3] * m.m[1][1] * m.m[2][0] + m.m[0][1] * m.m[1][0] * m.m[2][3] + m.m[0][0] * m.m[1][3] * m.m[2][1]) / determinant;

    num.m[3][0] = (-m.m[1][0] * m.m[2][1] * m.m[3][2] - m.m[1][1] * m.m[2][2] * m.m[3][0] - m.m[1][2] * m.m[2][0] * m.m[3][1] + m.m[1][2] * m.m[2][1] * m.m[3][0] + m.m[1][1] * m.m[2][0] * m.m[3][2] + m.m[1][0] * m.m[2][2] * m.m[3][1]) / determinant;
    num.m[3][1] = (m.m[0][0] * m.m[2][1] * m.m[3][2] + m.m[0][1] * m.m[2][2] * m.m[3][0] + m.m[0][2] * m.m[2][0] * m.m[3][1] - m.m[0][2] * m.m[2][1] * m.m[3][0] - m.m[0][1] * m.m[2][0] * m.m[3][2] - m.m[0][0] * m.m[2][2] * m.m[3][1]) / determinant;
    num.m[3][2] = (-m.m[0][0] * m.m[1][1] * m.m[3][2] - m.m[0][1] * m.m[1][2] * m.m[3][0] - m.m[0][2] * m.m[1][0] * m.m[3][1] + m.m[0][2] * m.m[1][1] * m.m[3][0] + m.m[0][1] * m.m[1][0] * m.m[3][2] + m.m[0][0] * m.m[1][2] * m.m[3][1]) / determinant;
    num.m[3][3] = (m.m[0][0] * m.m[1][1] * m.m[2][2] + m.m[0][1] * m.m[1][2] * m.m[2][0] + m.m[0][2] * m.m[1][0] * m.m[2][1] - m.m[0][2] * m.m[1][1] * m.m[2][0] - m.m[0][1] * m.m[1][0] * m.m[2][2] - m.m[0][0] * m.m[1][2] * m.m[2][1]) / determinant;

    return num;
}
Vector3 Normalize(const Vector3& v)
{
    float Normalize;
    Vector3 num;
    Normalize = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    num.x = v.x / Normalize;
    num.y = v.y / Normalize;
    num.z = v.z / Normalize;
    return num;
}
Vector3 Cross(const Vector3& v1, const Vector3& v2)
{
    Vector3 num;
    num = { v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x };
    return num;
}

Matrix4x4 MakeRotateXMatrix(float radian)
{
    Matrix4x4 num;
    num = { 1, 0, 0, 0,
        0, std::cos(radian), std::sin(radian), 0,
        0, std::sin(-radian), std::cos(radian), 0,
        0, 0, 0, 1 };
    return num;
}
Matrix4x4 MakeRotateYMatrix(float radian)
{
    Matrix4x4 num;
    num = { std::cos(radian), 0, std::sin(-radian), 0,
        0, 1, 0, 0,
        std::sin(radian), 0, std::cos(radian), 0,
        0, 0, 0, 1 };
    return num;
}
Matrix4x4 MakeRotateZMatrix(float radian)
{
    Matrix4x4 num;
    num = { std::cos(radian), std::sin(radian), 0, 0,
        std::sin(-radian), std::cos(radian), 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1 };
    return num;
}

Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate)
{
    Matrix4x4 rotateX = MakeRotateXMatrix(rotate.x);
    Matrix4x4 rotateY = MakeRotateYMatrix(rotate.y);
    Matrix4x4 rotateZ = MakeRotateZMatrix(rotate.z);
    Matrix4x4 rotateXYZ = Multiply(rotateX, Multiply(rotateY, rotateZ));

    Matrix4x4 num;
    num.m[0][0] = scale.x * rotateXYZ.m[0][0];
    num.m[0][1] = scale.x * rotateXYZ.m[0][1];
    num.m[0][2] = scale.x * rotateXYZ.m[0][2];
    num.m[0][3] = 0.0f * 0.0f * 0.0f * 0.0f;
    num.m[1][0] = scale.y * rotateXYZ.m[1][0];
    num.m[1][1] = scale.y * rotateXYZ.m[1][1];
    num.m[1][2] = scale.y * rotateXYZ.m[1][2];
    num.m[1][3] = 0.0f * 0.0f * 0.0f * 0.0f;
    num.m[2][0] = scale.z * rotateXYZ.m[2][0];
    num.m[2][1] = scale.z * rotateXYZ.m[2][1];
    num.m[2][2] = scale.z * rotateXYZ.m[2][2];
    num.m[2][3] = 0.0f * 0.0f * 0.0f * 0.0f;
    num.m[3][0] = translate.x;
    num.m[3][1] = translate.y;
    num.m[3][2] = translate.z;
    num.m[3][3] = 1.0f;
    return num;
}
Matrix4x4 MakePrespectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip)
{
    Matrix4x4 num;
    num = { (1 / aspectRatio) * (1 / tanf(fovY / 2)), 0, 0, 0, 0, (1 / tanf(fovY / 2)), 0, 0, 0, 0, farClip / (farClip - nearClip), 1, 0, 0, (-nearClip * farClip) / (farClip - nearClip) };
    return num;
}
Matrix4x4 MakeOrthographicMatrix(float left, float top, float right, float bottom, float nearClip, float farClip)
{
    Matrix4x4 num;
    num = { 2 / (right - left), 0, 0, 0, 0, 2 / (top - bottom), 0, 0, 0, 0, 1 / (farClip - nearClip), 0, (left + right) / (left - right),
        (top + bottom) / (bottom - top),
        nearClip / (nearClip - farClip), 1 };
    return num;
}
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth)
{
    Matrix4x4 num;
    num = { width / 2, 0, 0, 0, 0, -(height / 2), 0, 0, 0, 0, maxDepth - minDepth, 0, left + (width / 2), top + (height / 2), minDepth, 1 };
    return num;
}

// 3.座標変換
Vector3 TransForm(const Vector3& vector, const Matrix4x4& matrix)
{
    Vector3 result;
    result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + 1.0f * matrix.m[3][0];
    result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + 1.0f * matrix.m[3][1];
    result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + 1.0f * matrix.m[3][2];
    float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + 1.0f * matrix.m[3][3];
    assert(w != 0.0f);
    result.x /= w;
    result.y /= w;
    result.z /= w;

    return result;
}

void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix)
{
    const float kGridHalfWidth = 2.0f;
    const uint32_t kSubdivison = 10;
    const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivison);
    // 奥から手前に
    for (uint32_t xIndex = 0; xIndex <= kSubdivison; ++xIndex) {
        // 終点始点
        Vector3 numa = { -kGridHalfWidth + (kGridEvery * xIndex), 0.0f, kGridHalfWidth };
        Vector3 numb = { -kGridHalfWidth + (kGridEvery * xIndex), 0.0f, -kGridHalfWidth };
        // スクリーン座標
        Vector3 ndcnumA = TransForm(numa, viewProjectionMatrix);
        Vector3 screennumA = TransForm(ndcnumA, viewportMatrix);
        Vector3 ndcnumB = TransForm(numb, viewProjectionMatrix);
        Vector3 screennumB = TransForm(ndcnumB, viewportMatrix);

        if (xIndex != 5) {
            Novice::DrawLine((int)screennumA.x, (int)screennumA.y, (int)screennumB.x, (int)screennumB.y, 0xAAAAAAFF);
        } else {
            Novice::DrawLine((int)screennumA.x, (int)screennumA.y, (int)screennumB.x, (int)screennumB.y, 0x000000FF);
        }
    }

    for (uint32_t zIndex = 0; zIndex <= kSubdivison; ++zIndex) {
        // 終点始点
        Vector3 numa = { kGridHalfWidth, 0.0f, -kGridHalfWidth + (kGridEvery * zIndex) };
        Vector3 numb = { -kGridHalfWidth, 0.0f, -kGridHalfWidth + (kGridEvery * zIndex) };
        // スクリーン座標
        Vector3 ndcnumA = TransForm(numa, viewProjectionMatrix);
        Vector3 screennumA = TransForm(ndcnumA, viewportMatrix);
        Vector3 ndcnumB = TransForm(numb, viewProjectionMatrix);
        Vector3 screennumB = TransForm(ndcnumB, viewportMatrix);

        if (zIndex != 5) {
            Novice::DrawLine((int)screennumA.x, (int)screennumA.y, (int)screennumB.x, (int)screennumB.y, 0xAAAAAAFF);
        } else {
            Novice::DrawLine((int)screennumA.x, (int)screennumA.y, (int)screennumB.x, (int)screennumB.y, 0x000000FF);
        }
    }
}

Vector3 Prependicular(const Vector3& vector)
{
    if (vector.x != 0.0f || vector.y != 0.0f) {
        return { -vector.y, vector.x, 0.0f };
    }
    return { 0.0f, -vector.z, vector.y };
}

bool IsCollision(const Segment& segment, const Triangle& triangle, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix)
{
    Vector3 v01 = { triangle.vertices[1].x - triangle.vertices[0].x, triangle.vertices[1].y - triangle.vertices[0].y, triangle.vertices[1].z - triangle.vertices[0].z };
    Vector3 v12 = { triangle.vertices[2].x - triangle.vertices[1].x, triangle.vertices[2].y - triangle.vertices[1].y, triangle.vertices[2].z - triangle.vertices[1].z };
    Vector3 v02 = { triangle.vertices[0].x - triangle.vertices[2].x, triangle.vertices[0].y - triangle.vertices[2].y, triangle.vertices[0].z - triangle.vertices[2].z };

    Vector3 crossN = Cross(v01, v12);

    // 法線N
    crossN = Normalize(crossN);

    float d = Dot(triangle.vertices[0], crossN);

    Vector3 start = TransForm(TransForm(segment.origin, viewProjectionMatrix), viewportMatrix);
    Vector3 end = TransForm(TransForm(add(segment.origin, segment.diff), viewProjectionMatrix), viewportMatrix);

    if (Dot(segment.diff, crossN) == 0.0f) {
        Novice::DrawLine(int(start.x), int(start.y), int(end.x), int(end.y), WHITE);
        false;
    }

    float t = (d - Dot(segment.origin, crossN)) / Dot(crossN, segment.diff);

    if (t < 0.0f || t > 1.0f) {
        Novice::DrawLine(int(start.x), int(start.y), int(end.x), int(end.y), WHITE);
        return false;
    }

    Vector3 p = add(segment.origin, Multiply(t, segment.diff));

    Vector3 v1p = { p.x - triangle.vertices[1].x, p.y - triangle.vertices[1].y, p.z - triangle.vertices[1].z };
    Vector3 v2p = { p.x - triangle.vertices[2].x, p.y - triangle.vertices[2].y, p.z - triangle.vertices[2].z };
    Vector3 v0p = { p.x - triangle.vertices[0].x, p.y - triangle.vertices[0].y, p.z - triangle.vertices[0].z };

    // 確定底辺を組んだベクトルと頂点と衝突点pを結んだベクトルのクロス積を取る
    Vector3 cross01 = Cross(v01, v1p);
    Vector3 cross12 = Cross(v12, v2p);
    Vector3 cross20 = Cross(v02, v0p);

    // 全ての小三角形のクロス積と法線が同じ方向に向いていたら衝突
    if (Dot(cross01, crossN) >= 0.0f && Dot(cross12, crossN) >= 0.0f && Dot(cross20, crossN) >= 0.0f) {
        Novice::DrawLine(int(start.x), int(start.y), int(end.x), int(end.y), RED);
    } else {
        Novice::DrawLine(int(start.x), int(start.y), int(end.x), int(end.y), WHITE);
    }
    return false;
}

void DrawTriangle(const Triangle& trianglr, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
{

    Triangle wvptriangle;

    wvptriangle.vertices[0] = TransForm(TransForm(trianglr.vertices[0], viewProjectionMatrix), viewportMatrix);
    wvptriangle.vertices[1] = TransForm(TransForm(trianglr.vertices[1], viewProjectionMatrix), viewportMatrix);
    wvptriangle.vertices[2] = TransForm(TransForm(trianglr.vertices[2], viewProjectionMatrix), viewportMatrix);

    Novice::DrawLine(int(wvptriangle.vertices[0].x), int(wvptriangle.vertices[0].y), int(wvptriangle.vertices[1].x), int(wvptriangle.vertices[1].y), color);
    Novice::DrawLine(int(wvptriangle.vertices[1].x), int(wvptriangle.vertices[1].y), int(wvptriangle.vertices[2].x), int(wvptriangle.vertices[2].y), color);
    Novice::DrawLine(int(wvptriangle.vertices[2].x), int(wvptriangle.vertices[2].y), int(wvptriangle.vertices[0].x), int(wvptriangle.vertices[0].y), color);
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{

    // ライブラリの初期化
    Novice::Initialize(kWindowTitle, 1280, 720);

    // キー入力結果を受け取る箱
    char keys[256] = { 0 };
    char preKeys[256] = { 0 };

    Triangle triangle;
    triangle.vertices[0] = { -1.0f, 0.0f, 0.0f };
    triangle.vertices[1] = { 0.0f, 1.0f, 0.0f };
    triangle.vertices[2] = { 1.0f, 0.0f, 0.0f };

    // 今回使うやつ
    Segment segment { { -2.0f, -1.0f, 0.0f }, { 3.0f, 2.0f, 2.0f } };

    Vector3 cameraTranslate { 0.0f, 1.9f, -6.49f };
    Vector3 cameraRotate { 0.26f, 0.0f, 0.0f };

    float kWindowWidth = 1280.0f;
    float kWindowHeight = 720.0f;

    Matrix4x4 worldMatrix = MakeAffineMatrix({ 1.0f, 1.0f, 1.0f }, { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f });
    Matrix4x4 cameraMatrix = MakeAffineMatrix({ 1.0f, 1.0f, 1.0f }, cameraRotate, cameraTranslate);
    Matrix4x4 viewMatrix = Inverse(cameraMatrix);
    Matrix4x4 projectionMatrix = MakePrespectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
    Matrix4x4 worldViewProjectionMatrix = Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix));
    Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

    Vector3 start = TransForm(TransForm(segment.origin, worldViewProjectionMatrix), viewportMatrix);
    Vector3 end = TransForm(TransForm(add(segment.origin, segment.diff), worldViewProjectionMatrix), viewportMatrix);

    // ウィンドウの×ボタンが押されるまでループ
    while (Novice::ProcessMessage() == 0) {
        // フレームの開始
        Novice::BeginFrame();

        // キー入力を受け取る
        memcpy(preKeys, keys, 256);
        Novice::GetHitKeyStateAll(keys);

        ///
        /// ↓更新処理ここから
        ///

        ImGui::Begin("window");
        ImGui::DragFloat3("CameraTranslate", &cameraTranslate.x, 0.01f);
        ImGui::DragFloat3("CameraRotate", &cameraRotate.x, 0.01f);
        ImGui::DragFloat3("segment_origin", &segment.origin.x, 0.01f);
        ImGui::DragFloat3("segment_diff", &segment.diff.x, 0.01f);
        ImGui::DragFloat3("triangle.vertices[0]", &triangle.vertices[0].x);
        ImGui::DragFloat3("triangle.vertices[1]", &triangle.vertices[1].x);
        ImGui::DragFloat3("triangle.vertices[2]", &triangle.vertices[2].x);
        ImGui::End();

        worldMatrix = MakeAffineMatrix({ 1.0f, 1.0f, 1.0f }, { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f });
        cameraMatrix = MakeAffineMatrix({ 1.0f, 1.0f, 1.0f }, cameraRotate, cameraTranslate);
        viewMatrix = Inverse(cameraMatrix);
        projectionMatrix = MakePrespectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
        worldViewProjectionMatrix = Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix));
        viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

        IsCollision(segment, triangle, worldViewProjectionMatrix, viewportMatrix);

        ///
        /// ↑更新処理ここまで
        ///

        ///
        /// ↓描画処理ここから
        ///

        DrawGrid(worldViewProjectionMatrix, viewportMatrix);
        DrawTriangle(triangle, worldViewProjectionMatrix, viewportMatrix, 0xFFFFFFFF);

        ///
        /// ↑描画処理ここまで
        ///

        // フレームの終了
        Novice::EndFrame();

        // ESCキーが押されたらループを抜ける
        if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
            break;
        }
    }

    // ライブラリの終了
    Novice::Finalize();
    return 0;
}
