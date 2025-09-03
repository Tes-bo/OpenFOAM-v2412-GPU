/*=========================================================================

   Program: ParaView
   Module:    pqRecordEventsDialog.cxx

   Copyright (c) 2005-2008 Sandia Corporation, Kitware Inc.
   All rights reserved.

   ParaView is a free software; you can redistribute it and/or modify it
   under the terms of the ParaView license version 1.2.

   See License_v1.2.txt for the full ParaView license.
   A copy of this license can be obtained by contacting
   Kitware Inc.
   28 Corporate Drive
   Clifton Park, NY 12065
   USA

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

// Qt includes
#include <QDebug>
#include <QFile>
#include <QFileDialog>
#include <QMainWindow>
#include <QPushButton>
#include <QTextStream>
#include <QTimer>

// QtTesting includes
#include "pqEventComment.h"
#include "pqEventRecorder.h"
#include "pqEventTypes.h"
#include "pqRecordEventsDialog.h"
#include "pqTestUtility.h"

#include "ui_pqRecordEventsDialog.h"

//////////////////////////////////////////////////////////////////////////////////
// pqImplementation

struct pqRecordEventsDialog::pqImplementation
{
public:
  pqImplementation(pqEventRecorder* recorder, pqTestUtility* testUtility)
    : Recorder(recorder)
    , TestUtility(testUtility)
  {
  }

  ~pqImplementation()
  {
    this->Recorder = 0;
    this->TestUtility = 0;
  }

  Ui::pqRecordEventsDialog Ui;

  pqEventRecorder* Recorder;
  pqTestUtility* TestUtility;
};

///////////////////////////////////////////////////////////////////////////////////
// pqRecordEventsDialog

// ----------------------------------------------------------------------------
pqRecordEventsDialog::pqRecordEventsDialog(
  pqEventRecorder* recorder, pqTestUtility* testUtility, QWidget* Parent)
  : QDialog(Parent)
  , Implementation(new pqImplementation(recorder, testUtility))
{
  this->Implementation->Ui.setupUi(this);
  this->setWindowFlags(this->windowFlags().setFlag(Qt::WindowContextHelpButtonHint, false));
  this->setWindowFlags(this->windowFlags().setFlag(Qt::WindowStaysOnTopHint, true));

  this->ignoreObject(this);

  this->setWindowTitle(tr("Recording User Input"));
  this->setObjectName("");

  QObject::connect(this->Implementation->TestUtility->eventTranslator(),
    SIGNAL(recordEvent(QString, QString, QString, int)), this,
    SLOT(onEventRecorded(QString, QString, QString, int)));

  QObject::connect(
    this->Implementation->Ui.commentAddButton, SIGNAL(clicked()), this, SLOT(addComment()));

  QObject::connect(this->Implementation->Ui.checkButton, SIGNAL(toggled(bool)),
    this->Implementation->Recorder, SLOT(check(bool)));

  QObject::connect(this->Implementation->Ui.recordPauseButton, SIGNAL(toggled(bool)),
    this->Implementation->Recorder, SLOT(unpause(bool)));

  QObject::connect(this->Implementation->Ui.continuousFlush, SIGNAL(toggled(bool)),
    this->Implementation->Recorder, SLOT(setContinuousFlush(bool)));

  if (this->Implementation->TestUtility->supportsDashboardMode())
  {
    QObject::connect(this->Implementation->Ui.dashboardMode, SIGNAL(toggled(bool)), this,
      SLOT(onDashboardModeToggled(bool)));
  }
  else
  {
    this->Implementation->Ui.dashboardMode->setVisible(false);
  }

  QObject::connect(this->Implementation->Ui.recordInteractionTimings, SIGNAL(toggled(bool)),
    this->Implementation->Recorder, SLOT(setRecordInteractionTimings(bool)));

  QObject::connect(this->Implementation->Recorder, SIGNAL(started()), this, SLOT(updateUi()));
  QObject::connect(this->Implementation->Recorder, SIGNAL(stopped()), this, SLOT(updateUi()));
}

// ----------------------------------------------------------------------------
pqRecordEventsDialog::~pqRecordEventsDialog()
{
  delete Implementation;
}

// ----------------------------------------------------------------------------
void pqRecordEventsDialog::ignoreObject(QObject* object)
{
  this->Implementation->TestUtility->eventTranslator()->ignoreObject(object);
  Q_FOREACH (QObject* child, object->children())
  {
    this->ignoreObject(child);
  }
}

// ----------------------------------------------------------------------------
void pqRecordEventsDialog::done(int value)
{
  this->Implementation->TestUtility->stopRecords(value);
  QDialog::done(value);
}

// ----------------------------------------------------------------------------
void pqRecordEventsDialog::onEventRecorded(
  const QString& widget, const QString& command, const QString& argument, int /*eventType*/)
{
  if (!this->Implementation->Recorder->isRecording())
  {
    return;
  }

  this->Implementation->Ui.eventWidgetEdit->setText(widget);
  this->Implementation->Ui.eventCommandEdit->setText(command);
  this->Implementation->Ui.eventArgumentEdit->setText(argument);
  int newValue = this->Implementation->Ui.nbEvents->value() + 1;
  this->Implementation->Ui.nbEvents->display(newValue);
}

// ----------------------------------------------------------------------------
void pqRecordEventsDialog::addComment()
{
  if (this->Implementation->Ui.blockingCheckBox->isChecked())
  {
    this->Implementation->Recorder->translator()->eventComment()->recordCommentBlock(
      this->Implementation->Ui.commentTextEdit->toPlainText());
  }
  else
  {
    this->Implementation->Recorder->translator()->eventComment()->recordComment(
      this->Implementation->Ui.commentTextEdit->toPlainText());
  }
  this->Implementation->Ui.commentTextEdit->clear();
}

// ----------------------------------------------------------------------------
void pqRecordEventsDialog::updateUi()
{
  bool recording = this->Implementation->Recorder->isRecording();
  this->Implementation->Ui.recordPauseButton->setChecked(recording);
  this->Implementation->Ui.dashboardMode->setEnabled(recording);
}

// ----------------------------------------------------------------------------
void pqRecordEventsDialog::onDashboardModeToggled(bool toggle)
{
  this->Implementation->TestUtility->setDashboardMode(toggle);

  const QWidgetList& list = QApplication::topLevelWidgets();
  QMainWindow* mainWindow = nullptr;
  for (QWidget* widg : list)
  {
    mainWindow = qobject_cast<QMainWindow*>(widg);
    if (mainWindow)
    {
      break;
    }
  }
  this->Implementation->Recorder->translator()->recordDashboardModeToggle(mainWindow, toggle);
}
