#include <math.h>
#include <float.h>
#include <complex>
#include <iostream>
#include <algorithm>
#include <queue>
#include "Polynomials-main\Polynomial.cpp"



/*
    Имплементация метода решения уравнения четвертой степени — The fast quartic solver
    Информация о методе — The_fast_quartic_solver-Strobach-2010.pdf
    Работу выполнил — Федоров Алексей (https://github.com/Alexis777777)
*/



//Быстрый алгоритм наименьших квадратов для вычисления коэффициентов {γ0, δ0}
//для заданного набора{ α0, β0 }
template <typename fp_t>
unsigned int least_squares(fp_t a, fp_t b, fp_t c, fp_t d
    , fp_t alpha0, fp_t beta0, vector<fp_t>& coeff) {

    fp_t F1, F2;   
    F1 = fma(beta0, beta0, fma(alpha0, alpha0, static_cast<fp_t>(1.0L))); //F1 = 1 + pow(alpha0, 2) + pow(beta0, 2);   
    F2 = fma(alpha0, beta0, alpha0);        //F2 = alpha0 * (1 + beta0);

    fp_t c1, c2;
    c1 = fma(alpha0, b - beta0, -alpha0) + fma(beta0, c, a);     //c1 = a - alpha0 + alpha0 * (b - beta0) + beta0 * c;
    c2 = fma(alpha0, c, b) + fma(beta0, d, -beta0);     //c2 = b - beta0 + alpha0 * c + beta0 * d;

    fp_t L1, L2, L3;
    L1 = sqrt(F1);
    L3 = F2 / L1;
    L2 = sqrt(fma(-F2 / F1, F2, F1));       //L2 = sqrt(F1 - (F2 / F1) * F2);

    fp_t y1, y2;
    y1 = c1 / L1;
    y2 = fma(-y1, L3, c2) / L2;     //y2 = (c2 - y1 * L3) / L2;

    //Находим коэффициенты γ0, δ0
    fp_t delta0, gamma0;
    delta0 = y2 / L2;
    gamma0 = fma(-delta0, L3, y1) / L1;     //gamma0 = (y1 - delta0 * L3) / L1;

    coeff[0] = delta0;
    coeff[1] = gamma0;

    return 0;
}



//Итеративный алгоритм уточнения коэффициентов {α, β, γ, δ} из заданного начального
//приближения{ α0, β0, γ0, δ0 }.
template <typename fp_t> 
unsigned int backward_optimizer(fp_t a, fp_t b, fp_t c, fp_t d
    , fp_t alpha0, fp_t beta0, fp_t gamma0, fp_t delta0, unsigned maxIters, vector<fp_t>& coeff) {

    fp_t x1, x2, x3, x4;
    fp_t y1, y2, y3, y4;

    fp_t U23, U33, U44, L43;

    //Начальное приближение коэффициентов.
    fp_t alpha = alpha0, beta = beta0, gamma = gamma0, delta = delta0;

    //Инициализация величин ошибок на коэффициентах.
    fp_t e1 = a - alpha - gamma;
    fp_t e2 = b - fma(alpha, gamma, beta + delta); //b - beta - alpha * gamma - delta;
    fp_t e3 = fma(-alpha, delta, fma(-beta, gamma, c)); //c - beta * gamma - alpha * delta;
    fp_t e4 = fma(-beta, delta, d); //d - beta * delta;

    //Величина суммарной ошибки
    fp_t eps = 0;
    //Счетчик итераций
    unsigned iters = 0;

    //Очередь для хранения значений eps для 4-х предыдущих итераций
    queue <fp_t> q; q.push(0); q.push(0); q.push(0); q.push(0);
    fp_t tmp1, tmp2, tmp3, tmp4;

    for (int i = 0; i < maxIters; i++) {

        //LU-факторизация.
        U23 = alpha - gamma;
        U33 = -fma(gamma, U23, delta - beta);   //beta - delta - gamma * U23;

        L43 = -delta * U23 / U33;

        U44 = -fma(L43, U23, delta - beta); //beta - delta - L43 * U23;

        //Вычисление компонент {x1, x2, x3, x3} вспомогательного вектора x. Lx = e → x.
        x1 = e1;
        x2 = fma(-gamma, x1, e2);   //e2 - gamma * x1;
        x3 = fma(-gamma, x2, fma(-delta, x1, e3));  //e3 - delta * x1 - gamma * x2;
        x4 = fma(-L43, x3, fma(-delta, x2, e4));    //e4 - delta * x2 - L43 * x3;

        //Определения компонент {y1, y2, y3, y3} вектора обновления y. Uy = x → y.
        y4 = x4 / U44;
        y3 = fma(-U23, y4, x3) / U33; //(x3 - U23 * y4) / U33;
        y2 = fma(-U23, y3, x2 - y4);    //x2 - U23 * y3 - y4;
        y1 = x1 - y3;

        //Обновление коэффициентов.
        alpha = alpha + y1;
        beta = beta + y2;
        gamma = gamma + y3;
        delta = delta + y4;

        //Вычисление величин ошибок на коэффициентах.
        e1 = a - alpha - gamma;
        e2 = b - fma(alpha, gamma, beta + delta);  //b - beta - alpha * gamma - delta;
        e3 = fma(-alpha, delta, fma(-beta, gamma, c));  //c - beta * gamma - alpha * delta;
        e4 = fma(-beta, delta, d);  //d - beta * delta;

        iters = i;

        //Величина суммарной ошибки
        eps = abs(e1) + abs(e2) + abs(e3) + abs(e4);        

        tmp1 = q.front(); q.pop();
        tmp2 = q.front(); q.pop();
        tmp3 = q.front(); q.pop();
        tmp4 = q.front(); q.pop();
        
        //Условие сходимости алгоритма.
        if (eps < 1e-07 or isEqual(eps, tmp1) or isEqual(eps, tmp2) or isEqual(eps, tmp3) or isEqual(eps, tmp4)) break;
        //eps < 1e-07 подобрано экспериментально

        //Обновление элементов очереди.
        q.push(tmp2); q.push(tmp3); q.push(tmp4); q.push(eps);
    }

    coeff[0] = alpha;
    coeff[1] = beta;
    coeff[2] = gamma;
    coeff[3] = delta;
    coeff[4] = eps;

    return iters;
}


//Функция вычисления коэффициентов квадратных уравнений
template<typename fp_t>
unsigned int FQScoeffs(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& coeffs,
        int maxIters = 10000) {

    vector<complex<fp_t>> x(4);
    ferrari_complex<fp_t>(n, a, b, c, d, x); //Получаем приближение корней методом Феррари

    //Сортируем корни в порядке убывания
    complex<fp_t> tmp = 0;
    for (int i = 0; i < 4; i++)
        for (int j = (4 - 1); j >= (i + 1); j--)
            if (abs(x[j]) > abs(x[j - 1])) { 
                tmp = x[j]; x[j] = x[j - 1]; x[j - 1] = tmp;
            }
   
    //Вычисляем начальное приближение коэффициентов {α01, β01} и {α02, β02}
    fp_t alpha01, beta01, gamma01, delta01;
    fp_t alpha02, beta02, gamma02, delta02;

    alpha01 = -(x[0] + x[1]).real(); beta01 = (x[0] * x[1]).real();
    alpha02 = -(x[1] + x[2]).real(); beta02 = (x[1] * x[2]).real();

    //С помощью алгоритма наименьших квадратов для каждого из набора найдем коэффициенты { γ01, δ01 } и { γ02, δ02 } соответственно.
    vector<fp_t> coeff1(2);
    vector<fp_t> coeff2(2);

    least_squares(a, b, c, d, alpha01, beta01, coeff1);
    least_squares(a, b, c, d, alpha02, beta02, coeff2);

    delta01 = coeff1[0];
    gamma01 = coeff1[1];

    delta02 = coeff2[0];
    gamma02 = coeff2[1];

    //Найдем коэффициенты для каждого из двух начальных приближений
    vector<fp_t> res1(5);
    vector<fp_t> res2(5);

    unsigned itersCount1 = backward_optimizer(a, b, c, d, alpha01, beta01, gamma01, delta01, maxIters, res1);
    unsigned itersCount2 = backward_optimizer(a, b, c, d, alpha02, beta02, gamma02, delta02, maxIters, res2);

    fp_t eps1 = res1[4];
    fp_t eps2 = res2[4];

    //Определим какой набор коэффициентов правильный
    fp_t alpha, beta, gamma, delta;

    if (itersCount1 < itersCount2) {
        alpha = res1[0]; beta = res1[1]; gamma = res1[2]; delta = res1[3];
    }
    else if (itersCount1 > itersCount2) {
        alpha = res2[0]; beta = res2[1]; gamma = res2[2]; delta = res2[3];
    }
    else if (eps1 < eps2) {
        alpha = res1[0]; beta = res1[1]; gamma = res1[2]; delta = res1[3];
    }
    else {
        alpha = res2[0]; beta = res2[1]; gamma = res2[2]; delta = res2[3];
    }

    coeffs[0] = alpha;
    coeffs[1] = beta;
    coeffs[2] = gamma;
    coeffs[3] = delta;
}



template<typename fp_t>
unsigned int FastQuarticSolver(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots) {

    // Нормировка коэффициентов
    if (isZero(n) || isinf(a /= n))
        return solveCubic(a, b, c, d, roots);
    if (isinf(b /= n))
        return 0;
    if (isinf(c /= n))
        return 0;
    if (isinf(d /= n))
        return 0;

    vector<complex<fp_t>> roots1(2);
    vector<complex<fp_t>> roots2(2);
    
    unsigned numberOfRoots = 0;

    // Объявление констант
    static const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    static const fp_t ONE_6TH = static_cast<fp_t>(1.0L / 6.0L);

    //Case 1. Два двойных корня / Один четвертичный корень(4 - кратный корень).
    fp_t alpha_ = a * ONE_HALF;
    fp_t beta_ = fma(-alpha_, alpha_, b) * ONE_HALF;

    fp_t E1 = fma(static_cast<fp_t>(-2.0L) * alpha_, beta_, c);
    fp_t E2 = fma(-beta_, beta_, d);

    if (isZero(E1) and isZero(E2)) {
        solveQuadratic<fp_t>(static_cast<fp_t>(1.0L), alpha_, beta_, roots1);
        
        //Определим действительные корни
        if (isZero(roots1[0].imag())) {
            roots[numberOfRoots] = roots1[0].real();    
            roots[numberOfRoots+1] = roots1[0].real();
            numberOfRoots += 2;
        }
        if (isZero(roots1[1].imag())) {
            roots[numberOfRoots] = roots1[1].real();
            roots[numberOfRoots+1] = roots1[1].real();
            numberOfRoots += 2;
        }

        return numberOfRoots;
    }
            
    //Case 2. Один тройной корень и один простой корень.
    alpha_ = a * ONE_HALF; beta_ = b * ONE_6TH;
    vector<fp_t> roots_(2);
    solveQuadratic<fp_t>(static_cast<fp_t>(1.0L), alpha_, beta_, roots_);
    
    fp_t x1 = roots_[0], x2 = -fma(static_cast<fp_t>(3.0L), x1, a);
    E1 = fma(x1 * x1, fma(static_cast<fp_t>(3.0L), x2, x1), c); //c + x1 * x1 * (x1 + 3 * x2);
    E2 = fma(-x1 * x1, x1 * x2, d);

    if (isZero(E1) and isZero(E2)) {
        roots = { x1, x2, x2, x2 }; return numberOfRoots += 4;
    }

    x1 = roots_[1]; x2 = -fma(static_cast<fp_t>(3.0L), x1, a);
    E1 = fma(x1 * x1, fma(static_cast<fp_t>(3.0L), x2, x1), c); //c + x1 * x1 * (x1 + 3 * x2);
    E2 = fma(-x1 * x1, x1 * x2, d);

    if (isZero(E1) and isZero(E2)) {
        roots = { x1, x2, x2, x2 }; return numberOfRoots += 4;
    }

    //Case 3. Один двойной корень и два простых корня / 4 простых корня. 
    vector<fp_t> coeffs(4);
    FQScoeffs(n, a, b, c, d, coeffs);
    fp_t alpha = coeffs[0], beta = coeffs[1], gamma = coeffs[2], delta = coeffs[3];

    solveQuadratic<fp_t>(static_cast<fp_t>(1.0L), alpha, beta, roots1);
    solveQuadratic<fp_t>(static_cast<fp_t>(1.0L), gamma, delta, roots2);

    //Определим действительные корни
    if (isZero(roots1[0].imag())) {       
        roots[numberOfRoots] = roots1[0].real();    numberOfRoots++;
    }
    if (isZero(roots1[1].imag())) {
        roots[numberOfRoots] = roots1[1].real();    numberOfRoots++; 
    }
    if (isZero(roots2[0].imag())) {
        roots[numberOfRoots] = roots2[0].real();    numberOfRoots++; 
    }
    if (isZero(roots2[1].imag())) {
        roots[numberOfRoots] = roots2[1].real();    numberOfRoots++; 
    }

    return numberOfRoots;
}
