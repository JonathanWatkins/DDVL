////////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////////

// Standard includes
#include <iostream>

// Own includes
//#include "definitions.h"

// Qt includes
#include <QMessageBox>
#include <QCloseEvent>
#include <QTimer>

// Widgets
#include "mdmainwin.h"
#include "ui_mdmainwin.h"
#include "GLWidget.hpp"
#include "CSimulation.hpp"

////////////////////////////////////////////////////////////////
// CONSTRUCTOR & DESTRUCTOR
////////////////////////////////////////////////////////////////

mdmainwin::mdmainwin(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::mdmainwin)
{
    // Set up user interface
    ui->setupUi(this);

	  // Set up simulation output text browser
    int font_size = 7;
    int tab_width = 8;
    // Create font for the text browser
    QFont log_font("Courier New", font_size + 2, QFont::Normal, false); /* TODO: Why does Qt remove 2 from the font size?? */
    // Create text browser
    QTextEdit *tb = ui->simulation_output_tb;
    tb->setAutoFormatting(QTextEdit::AutoNone); /* Don't format inserted text */
    tb->setFont(log_font);
    tb->setLineWrapMode(QTextEdit::WidgetWidth); /* Wrap text at widget edge */
    tb->setOverwriteMode(false); /* Don't overwrite other text when inserting new one */
    tb->setReadOnly(true); /* Don't accept user inputed text */
    tb->setTabChangesFocus(true); /* Change focus when tab is pressed */
    tb->setTabStopWidth(tab_width * font_size); /* Tab width in pixels */
    tb->setUndoRedoEnabled(false); /* Don't allow undo */
    tb->setWordWrapMode(QTextOption::WrapAnywhere); /* Use all characters places on each line */

	  // Start simulation directly when application has finished loading
    QTimer::singleShot(0, this, SLOT(on_start_simulation_pb_clicked()));
}

mdmainwin::~mdmainwin()
{
    
    delete simGL;
    //delete simulation;
    delete ui;
}

////////////////////////////////////////////////////////////////
// PRIVATE SLOTS
////////////////////////////////////////////////////////////////

void mdmainwin::on_start_simulation_pb_clicked()
{
		if (simulation->is_running()) 
		{
        // Inform the user that an operation is currently going on
        QMessageBox msg_box;
        msg_box.setText("An operation is currently being executed.");
        msg_box.setInformativeText("Please wait until the current operation has finished.");
        msg_box.setStandardButtons(QMessageBox::Ok);
        msg_box.setDefaultButton(QMessageBox::Ok);
        msg_box.exec();
        return;
    }
    else
    {
			simulation = new CSimulation();
		
			//ui->verticalLayout_6->addWidget(simulation);
		
			//simulation->resize(600,400);
		
			//simulation->show();
		
			std::cout << "Start Simulation Clicked\n";
			ui->statusbar->showMessage("Start Simulation Clicked");
		}
}

void mdmainwin::closeEvent(QCloseEvent *event)
{
  
}

void mdmainwin::static_process_events(void* void_ptr_mainwin)
{
  
}
