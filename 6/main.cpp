#include <hprint>

template <typename T, std::size_t N>
void print_vector(const std::string_view vector_name, const std::array<T, N>& values)
{
    hp::print<'-'>(vector_name, 36);

    std::print("{:>8} {:>18}\n", "i", "val");

    for (std::size_t element_index = 0; element_index < N; ++element_index)
    {
        std::print("{:>8} {:>18.10f}\n", element_index, values[element_index]);
    }

    std::print("\n");
}

// Modyfikacja głównej przekątnej w miejscu -> 1 krok algorytmu Thomasa
template <typename T, std::size_t N>
void modify_diagonal_in_place(std::array<T, N>& main_diagonal, const std::array<T, N - 1>& upper_diagonal, const std::array<T, N - 1>& lower_diagonal)
{
    for (std::size_t row_index = 1; row_index < N; ++row_index)
    {
        main_diagonal[row_index] -= lower_diagonal[row_index - 1] * upper_diagonal[row_index - 1] / main_diagonal[row_index - 1];
    }
}

// Rozwiązanie układu na podstawie zmodyfikowanej przekątnej -> 2 krok algorytmu Thomasa
template <typename T, std::size_t N>
void solve_tridiagonal_system(const std::array<T, N - 1>& upper_diagonal, const std::array<T, N>& modified_main_diagonal, const std::array<T, N - 1>& lower_diagonal, std::array<T, N>& right_side)
{
    // Podstawianie do przodu
    for (std::size_t row_index = 1; row_index < N; ++row_index)
    {
        right_side[row_index] -= lower_diagonal[row_index - 1] * right_side[row_index - 1] / modified_main_diagonal[row_index - 1];
    }

    // Ostatni element rozwiązania manualnie
    right_side[N - 1] /= modified_main_diagonal[N - 1];

    // Podstawianie do tyłu
    for (std::size_t row_index = N - 1; row_index-- > 0;)
    {
        right_side[row_index] = ( right_side[row_index] - upper_diagonal[row_index] * right_side[row_index + 1] ) / modified_main_diagonal[row_index];
    }
}

int main()
{
    constexpr std::size_t system_size = 5;

    // Dla macierzy trójdiagonalnej przechowujemy tylko niezerowe przekątne
    constexpr std::array<long double, system_size - 1> upper_diagonal{ -1, -3, 5, -7 };
    std::array<long double, system_size>               main_diagonal{ 100, 200, 300, 200, 100 };
    constexpr std::array<long double, system_size - 1> lower_diagonal{ 2, 4, -6, -8 };
    std::array<long double, system_size>               right_side{ 199, 195, 929, 954, 360 };

    print_vector("Upper diagonal (u)", upper_diagonal);
    print_vector("Main diagonal (d)", main_diagonal);
    print_vector("Lower diagonal (l)", lower_diagonal);
    print_vector("Right-hand side (b)", right_side);

    modify_diagonal_in_place(main_diagonal, upper_diagonal, lower_diagonal);

    print_vector("Modified main diagonal (eta)", main_diagonal);

    solve_tridiagonal_system(upper_diagonal, main_diagonal, lower_diagonal, right_side);

    print_vector("Solution vector (x)", right_side);
}
