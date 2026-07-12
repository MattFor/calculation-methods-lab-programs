#include <print>
#include <cmath>

// Żeby wszystko było równo
constexpr int LABEL_WIDTH = 45;

/* IEEE754
 * float: 1.11 * 10^-7
 * 1 znaku
 * 8 wykładnika
 * 23 mantysy
 *
 * double: 1.11 * 10^-16
 * 1 znaku
 * 11 wykładnika
 * 52 mantysy
 *
 * long double: 1.11 * 10^-19 [teoretycznie]
 * 1 znaku
 * 15 wykładnika
 * 1 + 63 mantysy
*/

template <typename T>
void epsilon(const std::string& name)
{
    // Na później taki pretty print
    auto print_row = [&](std::string_view label, auto value)
    {
        std::print("{:<{}} = {}\n", label, LABEL_WIDTH, value);
    };

    std::print("{:-^30}\n", std::format(" {} ", name));

    // Zmienna reprezentująca liczbę 1 w typie T
    T one(1);

    // Start od eps = 1, które będzie zmieszane aż znaleziona będzie najmniejsza wartość, dla której 1 + eps jest różne od 1
    T eps(1);

    // Szukany jest "machine epsilon"
    // W każdej iteracji sprawdzane jest, czy 1 + eps/2 nadal daje wynik większy od 1,
    // Jeśli tak to znaczy, że eps jest jeszcze za duże, więc jest dalej zmniejszany x2
    // Pętla kończy się, gdy eps/2 jest już tak małe, że dodanie go do 1 nie zmienia wartości (ze względu na ograniczoną precyzję reprezentacji zmiennoprzecinkowej)
    while (one + eps / T(2) > one)
    {
        eps /= T(2);
    }

    // W tym momencie eps jest najmniejszą liczbą taką,
    // że 1 + eps > 1 -> czyli machine epsilon dla typu T

    // Epsilon wyliczony algorytmicznie
    std::print("Wyliczony epsilon\t\t= {:.{}e}\n", eps, std::numeric_limits<T>::max_digits10);

    // Epsilon podany przez bibliotekę standardową
    std::print("Epsilon z biblioteki\t= {:.{}e}\n", std::numeric_limits<T>::epsilon(), std::numeric_limits<T>::max_digits10);
    // ----------------------------
    // Wyliczenie liczby bitów mantysy
    // W IEEE754 zachodzi zależność "eps = 2^(1 - p)" gdzie p = liczba bitów precyzji (mantysa + bit ukryty)
    // po przekształceniu log2(eps) = 1 - p -> p = 1 - log2(eps)
    const long double log2_eps = std::log2(static_cast<long double>(eps));

    // Oszacowanie liczby bitów precyzji
    const long double p_est = 1. - log2_eps;

    // Zaokrąglenie do najbliższej liczby całkowitej
    const long double p_round = std::roundl(p_est);

    // Liczba bitów mantysy bez bitu ukrytego
    const long double mantissa_bits_est = p_round - 1;

    // Eps mnożony przez 10 aż wartość >= 1
    // Liczba mnożeń = najmniejsze k takie, że eps * 10^k >= 1,
    // Czyli k = ceil(-log10(eps))
    T   v_simple      = eps;
    int manual_digits = 0;

    while (v_simple < T(1) && manual_digits < 10000) // Guard na infloop
    {
        v_simple *= T(10); // Przecinek w prawo o 1
        ++manual_digits;   // Jedno przesunięcie, jedna cyfra
    }

    const long double dec_digits_est = -llround(std::log10(static_cast<long double>(eps)));

    // manual_digits to liczba mnożeń, czyli liczba znaczących cyfr dziesiętnych maxymalna
    print_row(std::format("numeric_limits<{}>::digits10", name), std::numeric_limits<T>::digits10);
    print_row(std::format("numeric_limits<{}>::max_digits10", name), std::numeric_limits<T>::max_digits10);
    print_row("Znaczące cyfry (ręczne)", manual_digits);

    auto       auto_digits       = static_cast<long long>(std::floor(dec_digits_est));
    const auto auto_digits2_ceil = static_cast<long long>(std::ceil(dec_digits_est));
    const auto auto_digits3_auto = static_cast<long long>(dec_digits_est);

    // Przybliżenie -> liczba cyfr ≈ -log10(eps), ponieważ epsilon określa minimalną względną różnicę między reprezentowanymi liczbami
    print_row("Znaczące cyfry (automatycznie)", std::format("{} [zgadza się? {}]", auto_digits, auto_digits == manual_digits ? "TAK" : "NIE"));
    print_row("Znaczące cyfry dokładnie (od prawej zaniżone)", std::format("{}", mantissa_bits_est - auto_digits));
    print_row("Znaczące cyfry dokładnie (od prawej zawyżone)", std::format("{}", mantissa_bits_est - auto_digits2_ceil));
    print_row("Znaczące cyfry dokładnie (od prawej auto)", std::format("{}", mantissa_bits_est - auto_digits3_auto));
    std::print("{:<{}} = {}\n", "Wyliczony numer bitów mantysy (ułamek)", LABEL_WIDTH, mantissa_bits_est);

    T value = eps;

    // Dobierany typ całkowity o takim samym rozmiarze jak typ zmiennoprzecinkowy
    // float -> uint32_t
    // double -> uint64_t
    // long double (16B) -> uint128
    using UInt = std::conditional_t<sizeof(T) == 4, std::uint32_t, std::conditional_t<sizeof(T) == 8, std::uint64_t, std::conditional_t<sizeof(T) == 16, unsigned __int128, void>>>;

    // Jeśli typ ma inny rozmiar, error i przerwanie
    static_assert(!std::is_same_v<UInt, void>, "Nieobsługiwany rozmiar zmiennoprzecinkowy!");

    // Liczba bitów w typie
    constexpr std::size_t bits = sizeof(T) * 8;

    // Reinterpretacja bitów liczby zmiennoprzecinkowej
    // jako liczbę całkowitą
    auto raw = std::bit_cast<UInt>(value);

    // Jeśli liczba mieści się w 64 bitach to można ją wydrukować bezpośrednio
    std::print("\nReprezentacja raw binarna:\n");
    if constexpr (sizeof(UInt) <= 8)
    {
        std::print("{:0{}b}\n", raw, bits);
    }
    else
    {
        // Dla 128-bitów liczba dzielona jest na górne i dolne 64 bity
        auto high = static_cast<std::uint64_t>(raw >> 64);
        auto low  = static_cast<std::uint64_t>(raw);
        std::print("{:064b}{:064b}\n", high, low);
    }

    std::print("\n\n");
}

int main()
{
    epsilon<float>("float");
    epsilon<double>("double");
    epsilon<long double>("long double");
}