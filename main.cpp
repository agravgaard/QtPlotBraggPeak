#include <QApplication>
#include <QtGui/QPainter>
#include <QtPrintSupport/QPrinter>

#include <algorithm> // transform, for_each, fill, generate
#include <array>     // array
#include <cmath>     // M_PI
#include <iostream>  // cout, cerr, endl, <<
#include <memory>    // unique_pointer
#include <numeric>   // accumulate
#include <string>    // string, <<
#include <valarray>  // valarray
#include <vector>    // vector

#include "paraboliccylinder.h"
#include "qcustomplot.h"

// bethe for water and water equivalent
double bethe(const double T) {
  // const double pi = 3.1415926535897932384626433832795028841971693993751;
  // const double N_A = 6.022140857e23; // mol^-1
  // const double mass_to_energy = 931.4940954; //MeV
  // const double Z = 0.555086707; // Z/A of water
  // const double I = 67.2e-6; // Exitation potential of water [MeV]

  // classical electron radius = e^2/(4*pi*eps * m_e)
  const double r = 2.8179403227e-15 * 100.0; // convert to cm
  const double m_e_x2 =
      0.5109989461 * 2.0; // MeV times two to shave off some instructions

  // Convert beam mass from AMU to MeV
  const double M_b = 1.0072766 * 931.4940954;

  // Beta^2
  const double b_2 = T * (T + 2.0 * M_b) / ((T + M_b) * (T + M_b));
  // Beta^2 * Gamma^2
  const double b_2_x_g_2 = b_2 / (1.0 - b_2);

  //                  pi
  const double first = 2.0 * 3.14159265358979323846 * (r * r) * m_e_x2;

  //                    N_A [mol^-1]     Z/A of water
  const double second = 6.022140857e23 * 0.555086707 / b_2; // / M_m=1;

  //                                            I [MeV]
  const double third = log(m_e_x2 * b_2_x_g_2 / 67.2e-6) - b_2;

  return first * second * third * 0.1; // MeV/mm
}

double BortfeldBraggPeak(const double R0, const double phi0,
                         const double epsilon, const double sig,
                         const double z) {
  // This is the very specialised equation for water only. It's equation 28/29
  // in the paper
  const double toGray = 1.602E-10;
  if (z < (R0 - 10 * sig)) {
    const auto fac = (phi0) / (1.0 + 0.012 * R0);
    const auto term1 = 17.93 * pow(R0 - z, -0.435);
    const auto term2 = ((0.444 + 31.7 * epsilon) / R0) * pow(R0 - z, 0.565);
    return fac * (term1 + term2) * toGray;
  } else if (z < (R0 + 5 * sig)) {
    const auto D565 = dv(-0.565, -((R0 - z) / sig));
    const auto D1565 = dv(-1.565, -((R0 - z) / sig));

    const auto frontfac =
        ((exp((-pow(R0 - z, 2)) / (4.0 * pow(sig, 2))) * pow(sig, 0.565)) /
         (1.0 + 0.012 * R0)) *
        phi0;
    const auto bracfac =
        11.26 * D565 / sig + ((0.157 + 11.26 * epsilon) / R0) * D1565;
    return frontfac * bracfac * toGray;
  }
  return 0.0;
}

bool generate_dEdx(std::valarray<double> &vec, const double init_energy,
                   const double step_length) {
  if (init_energy < 0.0 || step_length <= 0.0) {
    return false;
  }

  double energy = init_energy;
  std::generate(std::begin(vec), std::end(vec), [&energy, step_length]() {
    double dEdx = bethe(energy);
    energy -= dEdx * step_length;
    return dEdx;
  });
  return true;
}

bool generate_BraggPeak(std::valarray<double> &vec,
                        const std::valarray<double> &x_points,
                        const double init_energy, const double weight) {
  if (x_points.size() == 0 || init_energy < 0) {
    return false;
  }
  const auto phi = 4.0e9 * weight; // Fluence: particles/cm2
  const auto eps = 0.1;      // fraction of primary fluence contribution to tail
  const auto p = 1.77;       // exponent of range-energy relation
  const auto alpha = 0.0022; // proportionality factor in cm * MeV^-p

  const auto R0 = alpha * pow(init_energy, p); // approx range in cm
  const auto sigma_mono =
      0.012 * pow(R0, 0.935); // width of gaussian range straggling in cm
  const auto dRdE = alpha * p * pow(init_energy, p - 1.0);
  const auto sigma_E0 =
      0.01 * init_energy; // width of gaussian energy spectrum in MeV
  const auto sigma = sqrt(pow(sigma_mono, 2) + pow(sigma_E0 * dRdE, 2));

  std::transform(std::begin(x_points), std::end(x_points), std::begin(vec),
                 [R0, phi, eps, sigma](auto x_val) {
                   return BortfeldBraggPeak(R0, phi, eps, sigma, x_val / 10.0);
                 });

  return true;
}

int main(int argc, char *argv[]) {
  QApplication app(argc, argv);

  const auto n_points = 1000;
  const auto step_length = 400.0 / n_points; // 40 cm
  std::valarray<double> x_points(n_points);  // mm
  auto i = -step_length;
  std::generate(std::begin(x_points), std::end(x_points), [&i, step_length]() {
    i += step_length;
    return i;
  });

  i = 0.0;
  std::valarray<double> y_points_dEdx(n_points);
  const double init_energy = 200.0; // MeV
  if (!generate_dEdx(y_points_dEdx, init_energy, step_length)) {
    std::cerr << "Wrong arguments, dumbass!" << std::endl;
    return -1;
  };

  std::valarray<double> y_points_dose_200(n_points);
  auto weight = 1.0;
  if (!generate_BraggPeak(y_points_dose_200, x_points, init_energy, weight)) {
    std::cerr << "Did you not fill x_points?" << std::endl;
    return -2;
  }

  std::cout << "SOBP time!" << std::endl;
  // Let's make a SOBP of N bragg peaks:
  std::array<std::valarray<double>, 6> peak_array;

  std::array<double, 6> weights = {{0.92, 0.35, 0.25, 0.2, 0.16, 0.14}};
  // {1.0, 0.4, 0.25, 0.2, 0.15, 0.1}}; good for -=10 MeV
  const auto weight_sum = std::accumulate(weights.begin(), weights.end(), 0.0);
  std::transform(weights.begin(), weights.end(), weights.begin(),
                 [weight_sum](auto val) { return val / weight_sum; });

  auto init_e = init_energy;

  std::transform(weights.begin(), weights.end(), peak_array.begin(),
                 [&init_e, &x_points](auto weight) {
                   std::valarray<double> vec(x_points.size());
                   generate_BraggPeak(vec, x_points, init_e, weight);
                   init_e -= 5.0;
                   return vec;
                 });

  std::cout << "Individual peaks created" << std::endl;

  std::valarray<double> sobp(peak_array.at(0).size());
  std::fill(std::begin(sobp), std::end(sobp), 0.0);
  std::for_each(peak_array.begin(), peak_array.end(),
                [&sobp](auto vec) { sobp += vec; });
  auto sobp_max = sobp.max();
  sobp /= sobp_max;
  std::transform(peak_array.begin(), peak_array.end(), peak_array.begin(),
                 [sobp_max](auto vec) { return vec / sobp_max; });
  y_points_dose_200 /= y_points_dose_200.max();
  std::cout << "SOBP created" << std::endl;

  // "Copy" valarrays to QVectors (compiler should make it a memmove)
  QVector<double> Qx_points(x_points.size());
  std::copy(std::begin(x_points), std::end(x_points), Qx_points.begin());

  QVector<double> Qy_points_dose_200(y_points_dose_200.size());
  std::copy(std::begin(y_points_dose_200), std::end(y_points_dose_200),
            Qy_points_dose_200.begin());

  QVector<double> Qy_points_dEdx(y_points_dEdx.size());
  std::copy(std::begin(y_points_dEdx), std::end(y_points_dEdx),
            Qy_points_dEdx.begin());

  QVector<double> Qsobp(sobp.size());
  std::copy(std::begin(sobp), std::end(sobp), Qsobp.begin());

  std::cout << "copy complete (mostly)" << std::endl;

  auto customPlot = std::make_unique<QCustomPlot>();
  // create graph and assign data to it:
  int n_graph = 0;

  /// start bragg peak
  customPlot->addGraph();
  customPlot->graph(n_graph)->addData(Qx_points, Qy_points_dose_200);
  customPlot->graph(n_graph)->setPen(QPen(Qt::blue));
  customPlot->graph(n_graph)->setName("200 MeV Bragg Peak");

  customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft |
                                                                  Qt::AlignTop);
  customPlot->legend->setVisible(true);

  // give the axes some labels:
  customPlot->xAxis->setLabel("Depth [mm]");
  customPlot->yAxis->setLabel("Dose/max(Dose)");
  // set axes ranges, so we see all data:
  customPlot->xAxis->setRange(0, n_points * step_length);
  customPlot->yAxis->setRange(0, 1.1);


  /// start dE/dx
  n_graph++;
  customPlot->yAxis2->setVisible(true);
  customPlot->yAxis2->setLabel("dE/dx [keV/µm]");
  // MeV/mm to keV / µm is 1:1
  customPlot->yAxis2->setRange(0, y_points_dEdx.max());
  customPlot->addGraph(customPlot->xAxis, customPlot->yAxis2);
  customPlot->graph(n_graph)->addData(Qx_points, Qy_points_dEdx);
  customPlot->graph(n_graph)->setPen(QPen(Qt::red));
  customPlot->graph(n_graph)->setName("Stopping power (dE/dx)");

  customPlot->replot();
  customPlot->savePdf("out_0.pdf");

  /// start weighted peaks
  auto first_weighted_peak = n_graph + 1;
  // for_each is meant for side effects
  std::for_each(peak_array.begin(), peak_array.end(),
                [&customPlot, &n_graph, &Qx_points](auto vec) {
                  QVector<double> Qvec(vec.size());
                  std::copy(std::begin(vec), std::end(vec), Qvec.begin());
                  n_graph++;
                  customPlot->addGraph();
                  customPlot->graph(n_graph)->addData(Qx_points, Qvec);
                  customPlot->graph(n_graph)->setPen(QPen(Qt::gray));
                  customPlot->graph(n_graph)->setName("Weighted Bragg Peaks");
                });

  std::for_each(weights.begin(), std::prev(weights.end()),
                [&customPlot, first_weighted_peak](auto) {
                  customPlot->legend->removeItem(first_weighted_peak);
                });

  /// start SOBP
  n_graph++;
  customPlot->addGraph();
  customPlot->graph(n_graph)->addData(Qx_points, Qsobp);
  customPlot->graph(n_graph)->setPen(QPen(Qt::black));
  customPlot->graph(n_graph)->setName("SOBP");


  customPlot->replot();
  std::string out_pdf("out_");
  out_pdf += std::to_string(weights.size());
  out_pdf += ".pdf";
  customPlot->savePdf(out_pdf.c_str());

  return app.exec();
}
