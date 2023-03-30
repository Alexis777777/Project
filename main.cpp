#include <iostream>
#include "FQS.cpp"

template<typename fp_t>
void testQuarticPolynomial(int testCount, long double maxDistance)
{
    unsigned P = 4; // Степень исходного полинома
    fp_t low = -1, high = 1; // Интервал на котором заданы корни полинома
    fp_t absMaxError, relMaxError; // Абсолютная и относительная погрешность по итогам пройденного теста
    fp_t absMaxErrorTotal = -1, relMaxErrorTotal = -1; // Итоговая максимальная абсолютная и относительная погрешность по итогам всех тестов
    long double absErrorAvg = 0, relErrorAvg = 0; // Средняя абсолютная и относительная погрешность по итогам всех тестов
    unsigned numberOfFoundRoots; // Количество найденных корней
    unsigned cantFind = 0; // Счетчик количества ситуаций, когда методу не удалось найти корни (numberOfFoundRoots == 0)
    vector<fp_t> coefficients(P + 1); // Вектор коэффициентов полинома
    unsigned count = 0; // Счетчик количества ситуаций, когда относительная погрешность больше определенного числа (relMaxError > n)
    int countExcessRoots = 0;
    int countLostRoots = 0;

    for (size_t i = 0; i < testCount; ++i)
    {
        vector<fp_t> foundRoots(P);
        vector<fp_t> trueRoots(P);
        int excessRoots = 0;
        int lostRoots = 0;

        //generate_polynomial<fp_t>(P, 0, 0, 0, static_cast<fp_t>(maxDistance), low, high, trueRoots, coefficients);
        generate_polynomial<fp_t>(P, 0, P, 0, static_cast<fp_t>(maxDistance), low, high, trueRoots, coefficients);
        numberOfFoundRoots = FastQuarticSolver<fp_t>(coefficients[4], coefficients[3], coefficients[2], coefficients[1], coefficients[0], foundRoots);

        if (numberOfFoundRoots > 0)
        {
            compare_roots<fp_t>(numberOfFoundRoots, P, foundRoots, trueRoots, absMaxError, relMaxError, excessRoots, lostRoots);

            absMaxErrorTotal = absMaxError > absMaxErrorTotal ? absMaxError : absMaxErrorTotal;
            absErrorAvg += absMaxError;

            relMaxErrorTotal = relMaxError > relMaxErrorTotal ? relMaxError : relMaxErrorTotal;
            relErrorAvg += relMaxError;

            countExcessRoots += excessRoots;
            countLostRoots += lostRoots;

            count += relMaxError > 1 ? 1 : 0;
        }
        else
        {
            countLostRoots += 4;
            cantFind += 1;
        }
    }

    absErrorAvg /= (testCount - cantFind);
    relErrorAvg /= (testCount - cantFind);

    if (PRINT)
    {
        cout << "QUARTIC TEST RESULTS" << endl;
        cout << "========================================" << endl;
        cout << "Max distance: " << maxDistance << endl;
        cout << "Total count of tests: " << testCount << endl;
        cout << "Couldn't find roots: " << cantFind << " times " << endl;
        cout << "----------------------------------------" << endl;
        cout << "Average absolute error: " << absErrorAvg << endl;
        cout << "Total maximum absolute error: " << absMaxErrorTotal << endl;
        cout << "Average relative error: " << relErrorAvg << endl;
        cout << "Total maximum relative error: " << relMaxErrorTotal << endl;
        cout << "----------------------------------------" << endl;
        cout << "Total count of lost roots: " << countLostRoots << endl;
        cout << "Total count of excess roots: " << countExcessRoots << endl;
        cout << "----------------------------------------" << endl;
        cout << "relMaxError > 1: " << count << " times" << endl;
        cout << "========================================" << endl;
    }
}


int main()
{
   testQuarticPolynomial<fp_t>(1000000, 1e-5);
}