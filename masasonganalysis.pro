#-------------------------------------------------
#
# Project created by QtCreator 2018-11-15T11:34:07
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
QT += charts
QT += multimedia

TARGET = masasonganalysis
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11

SOURCES += \
        main.cpp \
        mainwindow.cpp \
    qwidgetwavchart.cpp \
    qchartviewwav.cpp \
    spuce/filters/butterworth_allpass.cpp \
    spuce/filters/butterworth_fir.cpp \
    spuce/filters/butterworth_iir.cpp \
    spuce/filters/calculate_decimator_taps.cpp \
    spuce/filters/chebyshev2_iir.cpp \
    spuce/filters/chebyshev_iir.cpp \
    spuce/filters/create_remez_lpfir.cpp \
    spuce/filters/design_fir.cpp \
    spuce/filters/design_iir.cpp \
    spuce/filters/design_window.cpp \
    spuce/filters/elliptic_allpass.cpp \
    spuce/filters/elliptic_iir.cpp \
    spuce/filters/farrow_upsampler.cpp \
    spuce/filters/find_roots.cpp \
    spuce/filters/fir_coeff.cpp \
    spuce/filters/fir_inv_dft.cpp \
    spuce/filters/gaussian_fir.cpp \
    spuce/filters/iir_coeff.cpp \
    spuce/filters/raised_cosine_imp.cpp \
    spuce/filters/remez_estimate.cpp \
    spuce/filters/remez_fir.cpp \
    spuce/filters/root_raised_cosine_imp.cpp \
    spuce/filters/shelf_allpass1.cpp \
    spuce/filters/sinc_fir.cpp \
    spuce/filters/sinc_helper.cpp \
    spuce/filters/transform_fir.cpp \
    spuce/filters/window.cpp

HEADERS += \
        mainwindow.h \
    qwidgetwavchart.h \
    qchartviewwav.h \
    spuce/dsp_classes/circ_buffer.h \
    spuce/dsp_classes/delay.h \
    spuce/dsp_functions/convolve.h \
    spuce/dsp_functions/fliplr.h \
    spuce/dsp_functions/partial_convolve.h \
    spuce/filters/allpass.h \
    spuce/filters/allpass_1.h \
    spuce/filters/allpass_2nd.h \
    spuce/filters/biquad.h \
    spuce/filters/butterworth_allpass.h \
    spuce/filters/butterworth_fir.h \
    spuce/filters/butterworth_iir.h \
    spuce/filters/calculate_decimator_taps.h \
    spuce/filters/cascaded_cic.h \
    spuce/filters/chebyshev2_iir.h \
    spuce/filters/chebyshev_iir.h \
    spuce/filters/cic.h \
    spuce/filters/create_remez_lpfir.h \
    spuce/filters/cutboost.h \
    spuce/filters/decimator.h \
    spuce/filters/design_fir.h \
    spuce/filters/design_iir.h \
    spuce/filters/design_window.h \
    spuce/filters/elliptic_allpass.h \
    spuce/filters/elliptic_iir.h \
    spuce/filters/farrow.h \
    spuce/filters/farrow_upsampler.h \
    spuce/filters/find_roots.h \
    spuce/filters/fir.h \
    spuce/filters/fir_adapt.h \
    spuce/filters/fir_coeff.h \
    spuce/filters/fir_decim.h \
    spuce/filters/fir_interp.h \
    spuce/filters/fir_inv_dft.h \
    spuce/filters/gaussian_fir.h \
    spuce/filters/iir.h \
    spuce/filters/iir_1st.h \
    spuce/filters/iir_allpass1_sections.h \
    spuce/filters/iir_allpass1_sections_variable_delay.h \
    spuce/filters/iir_allpass_variable_cascade.h \
    spuce/filters/iir_coeff.h \
    spuce/filters/iir_comb.h \
    spuce/filters/iir_df.h \
    spuce/filters/iir_shelf.h \
    spuce/filters/lagrange.h \
    spuce/filters/notch_allpass.h \
    spuce/filters/notch_comb.h \
    spuce/filters/notch_iir.h \
    spuce/filters/raised_cosine.h \
    spuce/filters/raised_cosine_imp.h \
    spuce/filters/remez_estimate.h \
    spuce/filters/remez_fir.h \
    spuce/filters/root_raised_cosine.h \
    spuce/filters/root_raised_cosine_imp.h \
    spuce/filters/running_average.h \
    spuce/filters/running_sum.h \
    spuce/filters/scic.h \
    spuce/filters/shelf_allpass1.h \
    spuce/filters/sinc_fir.h \
    spuce/filters/sinc_helper.h \
    spuce/filters/transform_fir.h \
    spuce/filters/window.h \
    spuce/base_type.h \
    spuce/complex_operators.h \
    spuce/mixed_type.h \
    spuce/typedefs.h

FORMS += \
        mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target


win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../fftw-3.3.5-dll32/ -llibfftw3-3
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../fftw-3.3.5-dll32/ -llibfftw3-3d
else:unix: LIBS += -L$$PWD/../../fftw-3.3.5-dll32/ -llibfftw3-3

INCLUDEPATH += $$PWD/../../fftw-3.3.5-dll32
DEPENDPATH += $$PWD/../../fftw-3.3.5-dll32]

DISTFILES += \
    spuce/Doxyfile \
    spuce/CMakeLists.txt \
    spuce/Doxyfile \
    spuce/CMakeLists.txt
