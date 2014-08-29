#ifndef MDMAINWIN_H
#define MDMAINWIN_H

//Standard includes
//#include <vector>
#include <string>
#include <iostream>

// Own includes
//#include "mdsystem.h"
//#include "settings.h"

#include "CSimulation.hpp"
#include "GLWidget.hpp"

// Qt includes
#include <QMainWindow>

namespace Ui {
    class mdmainwin;
}

class mdmainwin : public QMainWindow
{
    Q_OBJECT

private:
    //Private variables
    //settings system_settings;

public:
    // Constructor and destructor
    explicit mdmainwin(QWidget *parent = 0);
    ~mdmainwin();

private slots:

    void closeEvent(QCloseEvent *event);

    // Push buttons
   
    void on_start_simulation_pb_clicked();
    

private:
    // Private functions
    //void write_to_text_browser(std::string output);

    // Static private functions
    //static void static_write_to_text_browser(void* void_ptr_mainwin, string output);
    static void static_process_events       (void* void_ptr_mainwin               );

    // Private variables
    //mdsystem simulation;

private:
    Ui::mdmainwin *ui;
    GLWidget* simGL;
    CSimulation* simulation;

};

#endif // MDMAINWIN_H
