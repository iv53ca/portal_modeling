#include <iostream>
#include <vector>

#include <QApplication>
#include <QBoxLayout>
#include <QDebug>
#include <QMainWindow>
#include <QVector>
#include <QWidget>

// include all the plotting stuff
#include <QChart>
#include <QtCharts>
#include <QChartView>
#include <QValueAxis>
#include <QLineSeries>

// placed Eigen as third-party alongside project
#include <Eigen>


double a31(double Ka, double l1, double l2, double J){
    return(-Ka * l1 / (J * (l1 + l2)));
}
double a32(double Ka, double l1, double l2, double J){
    return(Ka * l1 / (J * (l1 + l2)));
}
double a33(double b1, double L, double M, double l1, double l2, double J){
    return(-b1 * L / (M * (l1 + l2))) - b1 * L * pow(l1,2) / (J * (l1 + l2));
}
double a34(double b2, double L, double M, double l1, double l2, double J){
    return(-b2 * L / (M * (l1 + l2))) + b2 * L * l1 * l2 / (J * (l1 + l2));
}
double a35(double kt1, double L, double M, double l1, double l2, double J){
    return(kt1 * L / (M * (l1 + l2))) + kt1 * L * pow(l1, 2) / (J * (l1 + l2));
}
double a36(double kt2, double L, double M, double l1, double l2, double J){
    return(kt2 * L / (M * (l1 + l2))) - kt2 * L * l1 * l2 / (J * (l1 + l2));
}
double a41(double Ka, double l1, double l2, double J){
    return(Ka * l2 / (J * (l1 + l2)));
}
double a42(double Ka, double l1, double l2, double J){
    return(-Ka * l2 / (J * (l1 + l2)));
}
double a43(double b1, double L, double M, double l1, double l2, double J){
    return(-b1 * L / (M * (l1 + l2))) + b1 * L * l1 * l2 / (J * (l1 + l2));
}
double a44(double b2, double L, double M, double l1, double l2, double J){
    return(-b2 * L / (M * (l1 + l2))) - b2 * L * pow(l2, 2) / (J * (l1 + l2));
}
double a45(double kt1, double L, double M, double l1, double l2, double J){
    return(kt1 * L / (M * (l1 + l2))) - kt1 * L * l1 * l2 / (J * (l1 + l2));
}
double a46(double kt2, double L, double M, double l1, double l2, double J){
    return(kt2 * L / (M * (l1 + l2))) + kt2 * L * pow(l2, 2) / (J * (l1 + l2));
}
double j(double M, double l1, double l2){
    return(M * (pow(l1, 2) - l1 * l2 + pow(l2, 2)) / 3);
}


void plotGraph(QWidget *parent, const std::vector<double>& y1, const std::vector<double>& y2, const std::vector<double>& time) {
    if (!parent) {
        qDebug() << "plotGraph: Provide non-null parent";
        return;
    }

    if (!parent->layout()) {
        qDebug() << "plotGraph: Provide parent which has layout";
        return;
    }

    // Check that the input vectors are not empty
    if (y1.empty()  || y2.empty()  ||  time.empty()) {
        qDebug() << "Error: One or more input vectors are empty";
        return;
        }

    // Check that the input vectors have the same size
    if (y1.size() != time.size(), y2.size() != time.size(), y1.size() != y2.size()) {
        qDebug() << "Error: The input vectors have different sizes";
        return;
        }

    // Create a chart
    QChart *chart = new QChart();

    // Create x-axis (time)
    QValueAxis *axisX = new QValueAxis();

    // Create y-axis
    QValueAxis *axisY = new QValueAxis();

    // Set axis titles
    axisX->setTitleText("Time");
    axisY->setTitleText("y1");

    // Add axes to the chart
    chart->addAxis(axisX, Qt::AlignBottom);
    chart->addAxis(axisY, Qt::AlignLeft);

    // Create a line series
    QLineSeries *series1 = new QLineSeries();
    QLineSeries *series2 = new QLineSeries();

    // Fill the series with data from vectors
    for (size_t i = 0; i < y1.size(); ++i) {
        series1->append(time[i], y1[i]);
        series2->append(time[i], y2[i]);
    }

    // Add series to the chart
    chart->addSeries(series1);
    chart->addSeries(series2);

    // Attach series to axes
    series1->attachAxis(axisX);
    series1->attachAxis(axisY);
    series2->attachAxis(axisX);
    series2->attachAxis(axisY);

    // Create a chart view
    auto view = new QChartView(parent);
    view->setMinimumSize(600, 500);
    view->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);

    view->setChart(chart);

    parent->layout()->addWidget(view);
    view->show();
}


int main(int argc, char *argv[]) {
    //Parameters definition
    double midrule_crossbeam_mass = 25.0;
    double load_mass = 10.0;
    double crossbeam_length = 0.8;
    double damping_y1 = 5.0;
    double damping_y2 = 5.0;
    double stiffness = 52.52;
    double thrust_constant_y1 = 61.0;
    double thrust_constant_y2 = 61.0;
    double back_EMF_y1 = 49.6;
    double back_EMF_y2 = 49.6;
    double inductance_y1 = 5.07 / 1000;
    double inductance_y2 = 5.07 / 1000;
    double resistance_y1 = 8.4;
    double resistance_y2 = 8.4;

    //Example of different l1 + l2 = L
    double l1 = 0.3;
    double l2 = 0.5;

    //Exact parameters definition for current problem
    double J = j(midrule_crossbeam_mass, l1, l2);
    double A31 = a31(stiffness, l1, l2, J);
    double A32 = a32(stiffness, l1, l2, J);
    double A33 = a33(damping_y1, crossbeam_length, midrule_crossbeam_mass, l1, l2, J);
    double A34 = a34(damping_y2, crossbeam_length, midrule_crossbeam_mass, l1, l2, J);
    double A35 = a35(thrust_constant_y1, crossbeam_length, midrule_crossbeam_mass, l1, l2, J);
    double A36 = a36(thrust_constant_y2, crossbeam_length, midrule_crossbeam_mass, l1, l2, J);
    double A41 = a41(stiffness, l1, l2, J);
    double A42 = a42(stiffness, l1, l2, J);
    double A43 = a43(damping_y1, crossbeam_length, midrule_crossbeam_mass, l1, l2, J);
    double A44 = a44(damping_y2, crossbeam_length, midrule_crossbeam_mass, l1, l2, J);
    double A45 = a45(thrust_constant_y1, crossbeam_length, midrule_crossbeam_mass, l1, l2, J);
    double A46 = a46(thrust_constant_y2, crossbeam_length, midrule_crossbeam_mass, l1, l2, J);

    //Definition of matrices solution
    Eigen::Matrix <double, 6, 6> A{
            {0, 0, 1, 0, 0, 0},
            {0, 0, 0, 1, 0, 0},
            {A31, A32, A33, A34, A35, A36},
            {A41, A42, A43, A44, A45, A46},
            {0, 0, -back_EMF_y1 / inductance_y1, 0, -resistance_y1 / inductance_y1, 0},
            {0, 0, 0, -back_EMF_y2 / inductance_y2, 0, -resistance_y2 / inductance_y2}
    };

    Eigen::Matrix <double, 6, 2> B{
            {0, 0},
            {0, 0},
            {0, 0},
            {0, 0},
            {1 / inductance_y1, 0},
            {0, 1 / inductance_y2}
    };

    Eigen::Matrix <double, 2, 6> C{
            {1, 0, 0, 0, 0, 0},
            {0, 1, 0, 0, 0, 0}
    };

    Eigen::Matrix<double, 6, 1> x{
            {0},
            {0},
            {0},
            {0},
            {0},
            {0}
    };

    Eigen::Matrix<double, 2, 1> u{
            {1},
            {1}
    };


    double dt = 0.0001;
    double T = 1.0;

    int order = 4;
    Eigen::Matrix<double, 6, 1> x_prev[order];
    for (int i = 0; i < order; ++i) {
        x_prev[i] = x;
    }

    std::vector<double> y1(T / dt, 0);
    std::vector<double> y2(T / dt, 0);
    std::vector<double> time(T / dt, 0);
    int cnt = 0;

    //Linear multistep method of Adams-Bashforth-Moulton
    for (double t = 0; t < T; t += dt) {
        Eigen::Matrix<double, 6, 1> f_prev[order];
        for (int i = 0; i < order; ++i) {
            f_prev[i] = (A * x_prev[i] + B * u);
        }

        Eigen::Matrix<double, 6, 1> f = (A * x + B * u);

        x = x_prev[3] + dt / 24.0 * (55*f - 59*f_prev[0] + 37*f_prev[1] - 9*f_prev[2]);

        for (int i = 0; i < order - 1; ++i) {
            x_prev[i] = x_prev[i+1];
        }
        x_prev[order-1] = x;

        Eigen::Matrix<double, 2, 1> y = C * x;

        y1[cnt] = y[0];
        y2[cnt] = y[1];
        time[cnt] = t;

        cnt+=1;
    }

    QApplication app(argc, argv);

    QMainWindow mainWindow;
    mainWindow.setMinimumSize(800, 600);

    auto centralWidget = new QWidget(&mainWindow);
    mainWindow.setCentralWidget(centralWidget);

    auto layout = new QBoxLayout(QBoxLayout::TopToBottom, mainWindow.centralWidget());
    mainWindow.centralWidget()->setLayout(layout);

    mainWindow.show();

    plotGraph(mainWindow.centralWidget(), y1, y2,  time);

    return app.exec();
}
