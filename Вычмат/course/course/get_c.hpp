#pragma once
#include <vector>

using Matrix = std::vector<double>;

Matrix get_c();
Matrix create_inverse_matrix(Matrix m);
Matrix create_identity_matrix();
