/*=========================================================================

   Program: ParaView
   Module:    pqEventTranslator.cxx

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

#include "pqEventTranslator.h"

#include "pq3DViewEventTranslator.h"
#include "pqAbstractButtonEventTranslator.h"
#include "pqAbstractItemViewEventTranslator.h"
#include "pqAbstractSliderEventTranslator.h"
#include "pqBasicWidgetEventTranslator.h"
#include "pqCheckEventOverlay.h"
#include "pqComboBoxEventTranslator.h"
#include "pqDoubleSpinBoxEventTranslator.h"
#include "pqEventComment.h"
#include "pqEventTypes.h"
#include "pqLineEditEventTranslator.h"
#include "pqListViewEventTranslator.h"
#include "pqMenuEventTranslator.h"
#include "pqNativeFileDialogEventTranslator.h"
#include "pqObjectNaming.h"
#include "pqSpinBoxEventTranslator.h"
#include "pqTabBarEventTranslator.h"
#include "pqTableViewEventTranslator.h"
#include "pqTreeViewEventTranslator.h"

#include <QApplication>
#include <QCoreApplication>
#include <QElapsedTimer>
#include <QSet>
#include <QToolBar>
#include <QVector>
#include <QtDebug>

////////////////////////////////////////////////////////////////////////////////
// pqEventTranslator::pqImplementation
struct pqEventTranslator::pqImplementation
{
  pqImplementation()
  {
    this->EventComment = 0;
    this->Checking = false;
    this->Recording = false;

    // Create overlay and hide it
    this->CheckOverlay = new pqCheckEventOverlay(NULL);
    this->hideOverlay();

    this->RecordInteractionTimings = false;
  }

  ~pqImplementation()
  {
    if (this->EventComment)
    {
      delete this->EventComment;
    }
    delete this->CheckOverlay;
    this->CheckOverlayWidgetOn = NULL;
  }

  void hideOverlay()
  {
    // Hide the overlay
    this->CheckOverlay->hide();

    // Nullfied it's parent
    this->CheckOverlay->setParent(NULL);

    // Set the overlayed widget to null
    this->CheckOverlayWidgetOn = NULL;
  }

  pqEventComment* EventComment;
  /// Stores the working set of widget translators
  QList<pqWidgetEventTranslator*> Translators;
  /// Stores the set of objects that should be ignored when translating events
  QMap<QObject*, QRegularExpression> IgnoredObjects;

  // list of widgets for which mouse propagation will happen
  // we'll only translate the first and ignore the rest
  QList<QWidget*> MouseParents;

  // Checking Flag
  bool Checking;

  // Recording flag
  bool Recording;

  // Pointer to the overlay
  QPointer<pqCheckEventOverlay> CheckOverlay;

  // Pointer to the overlayed widget
  QPointer<QWidget> CheckOverlayWidgetOn;

  // Record interaction timings flag
  bool RecordInteractionTimings;

  // Timer to track time between user interactions
  QElapsedTimer InteractionsTimer;
};

////////////////////////////////////////////////////////////////////////////////
// pqEventTranslator

// ----------------------------------------------------------------------------
pqEventTranslator::pqEventTranslator(QObject* p)
  : QObject(p)
  , Implementation(new pqImplementation())
{
  // Ignore the overlay so it is transparent to events.
  this->ignoreObject(this->Implementation->CheckOverlay);
}

// ----------------------------------------------------------------------------
pqEventTranslator::~pqEventTranslator()
{
  this->stop();
  delete Implementation;
}

// ----------------------------------------------------------------------------
void pqEventTranslator::start()
{
  QCoreApplication::instance()->installEventFilter(this);
  Q_EMIT this->started();
}

// ----------------------------------------------------------------------------
void pqEventTranslator::stop()
{
  QCoreApplication::instance()->removeEventFilter(this);
  this->check(false);
  Q_EMIT this->stopped();
}

// ----------------------------------------------------------------------------
void pqEventTranslator::addDefaultWidgetEventTranslators(pqTestUtility* util)
{
  // Add generalistic translator first, then specific, in order for this to work
  addWidgetEventTranslator(new pqBasicWidgetEventTranslator());
  addWidgetEventTranslator(new pqAbstractButtonEventTranslator());
  addWidgetEventTranslator(new pqAbstractItemViewEventTranslator());
  addWidgetEventTranslator(new pqAbstractSliderEventTranslator());
  addWidgetEventTranslator(new pqComboBoxEventTranslator());
  addWidgetEventTranslator(new pqDoubleSpinBoxEventTranslator());
  addWidgetEventTranslator(new pqLineEditEventTranslator());
  addWidgetEventTranslator(new pqMenuEventTranslator());
  addWidgetEventTranslator(new pqSpinBoxEventTranslator());
  addWidgetEventTranslator(new pqTabBarEventTranslator());
  addWidgetEventTranslator(new pqTreeViewEventTranslator());
  addWidgetEventTranslator(new pqTableViewEventTranslator());
  addWidgetEventTranslator(new pqListViewEventTranslator());
  addWidgetEventTranslator(new pq3DViewEventTranslator("QGLWidget"));
  addWidgetEventTranslator(new pqNativeFileDialogEventTranslator(util));
}

// ----------------------------------------------------------------------------
void pqEventTranslator::addWidgetEventTranslator(pqWidgetEventTranslator* Translator)
{
  if (Translator)
  {
    // We Check if the translator has already been added previously.
    int index = this->getWidgetEventTranslatorIndex(Translator->metaObject()->className());
    if (index != -1)
    {
      return;
    }

    this->Implementation->Translators.push_front(Translator);
    Translator->setParent(this);

    // Legacy Connection, for translator not using event types
    QObject::connect(Translator, SIGNAL(recordEvent(QObject*, const QString&, const QString&)),
      this, SLOT(onRecordEvent(QObject*, const QString&, const QString&)));

    // Connect record event
    QObject::connect(Translator, SIGNAL(recordEvent(QObject*, const QString&, const QString&, int)),
      this, SLOT(onRecordEvent(QObject*, const QString&, const QString&, int)));

    // Connect resize specific overlay
    QObject::connect(Translator, SIGNAL(specificOverlay(const QRect&)), this,
      SLOT(setOverlayGeometry(const QRect&)));
  }
}

// ----------------------------------------------------------------------------
bool pqEventTranslator::removeWidgetEventTranslator(const QString& className)
{
  int index = this->getWidgetEventTranslatorIndex(className);
  if (index == -1)
  {
    return false;
  }

  this->Implementation->Translators.removeAt(index);
  return true;
}

// ----------------------------------------------------------------------------
pqWidgetEventTranslator* pqEventTranslator::getWidgetEventTranslator(const QString& className)
{
  int index = this->getWidgetEventTranslatorIndex(className);
  if (index == -1)
  {
    return 0;
  }

  return this->Implementation->Translators.at(index);
}

// ----------------------------------------------------------------------------
int pqEventTranslator::getWidgetEventTranslatorIndex(const QString& className)
{
  for (int i = 0; i < this->Implementation->Translators.count(); ++i)
  {
    if (this->Implementation->Translators.at(i)->metaObject()->className() == className)
    {
      return i;
    }
  }
  return -1;
}

// ----------------------------------------------------------------------------
QList<pqWidgetEventTranslator*> pqEventTranslator::translators() const
{
  return this->Implementation->Translators;
}

// ----------------------------------------------------------------------------
void pqEventTranslator::addDefaultEventManagers(pqTestUtility* util)
{
  this->Implementation->EventComment = new pqEventComment(util);
  QObject::connect(this->Implementation->EventComment,
    SIGNAL(recordComment(QObject*, QString, QString)), this,
    SLOT(onRecordEvent(QObject*, QString, QString)));
}

// ----------------------------------------------------------------------------
pqEventComment* pqEventTranslator::eventComment() const
{
  return this->Implementation->EventComment;
}

// ----------------------------------------------------------------------------
void pqEventTranslator::ignoreObject(QObject* object, QRegularExpression commandFilter)
{
  this->Implementation->IgnoredObjects.insert(object, commandFilter);
}

// ----------------------------------------------------------------------------
bool pqEventTranslator::eventFilter(QObject* object, QEvent* event)
{
  if (this->Implementation->Recording)
  {
    if (object->isWindowType())
    {
      return false;
    }

    // Only widgets
    QWidget* widget = qobject_cast<QWidget*>(object);
    if (widget != NULL)
    {

      // mouse events are propagated to parents
      // our event translators/players don't quite like that,
      // so lets consume those extra ones
      if (event->type() == QEvent::MouseButtonPress ||
        event->type() == QEvent::MouseButtonDblClick || event->type() == QEvent::MouseMove ||
        event->type() == QEvent::Enter || event->type() == QEvent::Leave ||
        event->type() == QEvent::MouseButtonRelease || event->type() == QEvent::ContextMenu)
      {
        if (!this->Implementation->MouseParents.empty() &&
          this->Implementation->MouseParents.first() == object)
        {
          // right on track
          this->Implementation->MouseParents.removeFirst();
          return false;
        }

        // find the chain of parent that will get this mouse event
        this->Implementation->MouseParents.clear();
        for (QWidget* w = widget->parentWidget(); w; w = w->parentWidget())
        {
          this->Implementation->MouseParents.append(w);
          if (w->isWindow() || w->testAttribute(Qt::WA_NoMousePropagation))
          {
            break;
          }
        }
      }

      // Checking mode
      if (this->Implementation->Checking)
      {
        // In Gl Case, parentless widget is not transparent to mouse event
        // The event is  applied to the overlayed widget or an another top widget
        // (before ignoredObjects)
        // TODO : use mask instead

        // Look behind the overlay widget, this is relevant only if :
        // There is already an overlay drawn
        // This is the gl case
        // This is a mouse event we manage
        // The overlay receive the event
        if (this->Implementation->CheckOverlayWidgetOn != NULL &&
          this->Implementation->CheckOverlay->GlWidget &&
          widget == this->Implementation->CheckOverlay &&
          (event->type() == QEvent::MouseButtonRelease || event->type() == QEvent::MouseMove))
        {

          // We are about to look for another top widget window behind the mouse cursor
          bool foundTop = false;
          QWidget* topWidget;
          // recover all top widgets
          QWidgetList topWidgets = QApplication::topLevelWidgets();
          Q_FOREACH (topWidget, topWidgets)
          {
            // only the visible ones
            if (!topWidget->isHidden())
            {
              // Check it is not the overlay, and it contains the mouse cursor
              if (topWidget != this->Implementation->CheckOverlay &&
                topWidget->geometry().contains(static_cast<QMouseEvent*>(event)->globalPos(), true))
              {
                // Recover the child widget onder the cursor, if any
                QWidget* childWidget = topWidget->childAt(
                  topWidget->mapFromGlobal(static_cast<QMouseEvent*>(event)->globalPos()));

                // If child exist, check it is not the overlayed widget and indeed a new widget
                if (childWidget == NULL ||
                  (childWidget != NULL &&
                      childWidget != this->Implementation->CheckOverlayWidgetOn))
                {
                  // if a child exist, use it
                  if (childWidget != NULL)
                  {
                    topWidget = childWidget;
                  }

                  // Position top widget flag
                  foundTop = true;
                  break;
                }
              }
            }
          }
          if (foundTop)
          {
            // If we found a top widget behin the cursor, use it
            widget = topWidget;
          }
          else
          {
            // If not, use the widget behind the overlay
            widget = this->Implementation->CheckOverlayWidgetOn;
          }
        }

        // Leaving the checkOverlay when a widget without grand parent
        // then hiding overlay (before ignoredObjects)
        if (this->Implementation->CheckOverlay == widget &&
          !this->Implementation->CheckOverlay->Specific)
        {
          if (event->type() == QEvent::Leave)
          {
            this->Implementation->hideOverlay();
          }
        }

        // Ignore object if specified
        if (this->Implementation->IgnoredObjects.contains(widget))
        {
          return false;
        }

        // Do not check QToolBar as they can generate buggy overlay
        if (qobject_cast<QToolBar*>(widget) != NULL)
        {
          return false;
        }

        // Mouse Move on a non-previously overlayed widget
        if (event->type() == QEvent::MouseMove &&
          this->Implementation->CheckOverlayWidgetOn != widget)
        {
          // Check for any valid translator
          bool validTranslator = false;
          bool error;
          for (int i = 0; i != this->Implementation->Translators.size(); ++i)
          {
            // This will not record a check event, only check if a translator can record a check
            // event
            if (this->Implementation->Translators[i]->translateEvent(
                  widget, event, pqEventTypes::CHECK_EVENT, error))
            {
              validTranslator = true;
              break;
            }
          }

          // Draw overlay
          this->Implementation->hideOverlay();

          // Only on non-window widget
          // TODO handle the overlay on window widget
          if (!widget->isWindow() && widget->parent() != NULL)
          {
            // Check if widget is glWidget
            this->Implementation->CheckOverlay->GlWidget = widget->inherits("QVTKWidget");

            // Check if widget is parent to gl widget
            QList<QWidget*> children = widget->findChildren<QWidget*>();
            Q_FOREACH (QWidget* child, children)
            {
              if (child->inherits("QVTKWidget"))
              {
                // Ignore widget containing a Gl widget as a child
                this->ignoreObject(widget);
                return false;
              }
            }

            // Check widget size is valid, if not do not draw overlay
            // Widget can still be clicked and checked
            int sizeThreshold =
              pqCheckEventOverlay::OVERLAY_MARGIN + 2 * pqCheckEventOverlay::OVERLAY_PEN_WIDTH;
            if (widget->width() < sizeThreshold || widget->height() < sizeThreshold)
            {
              return false;
            }

            // Set the validity of the overlay
            this->Implementation->CheckOverlay->Valid = validTranslator;

            // Set parent of the overlay to be parent of the overlayed widget
            this->Implementation->CheckOverlay->setParent(qobject_cast<QWidget*>(widget->parent()));

            if (this->Implementation->CheckOverlay->GlWidget)
            {
              // Cannot draw QPainter directive in openGl context, bust use another context, aka
              // another window
              // this->Implementation->CheckOverlay->setWindowFlags(Qt::ToolTip |
              // Qt::FramelessWindowHint); // ToolTip is always on top
              this->Implementation->CheckOverlay->setWindowFlags(
                Qt::Tool | Qt::FramelessWindowHint); // Tool generate avatar
              this->Implementation->CheckOverlay->setAttribute(Qt::WA_NoSystemBackground, true);
              this->Implementation->CheckOverlay->setAttribute(Qt::WA_TranslucentBackground, true);
              this->Implementation->CheckOverlay->setAttribute(Qt::WA_PaintOnScreen, true);

              // Resize and move widget
              QRect geometry = widget->geometry();
              geometry.moveTo(widget->mapToGlobal(QPoint(0, 0)));
              this->setOverlayGeometry(geometry, false);
            }
            else
            {
              // Restore window flags
              this->Implementation->CheckOverlay->setWindowFlags(Qt::Widget);
              this->Implementation->CheckOverlay->setAttribute(Qt::WA_NoSystemBackground, false);
              this->Implementation->CheckOverlay->setAttribute(Qt::WA_TranslucentBackground, false);
              this->Implementation->CheckOverlay->setAttribute(Qt::WA_PaintOnScreen, false);

              // Set overlay geometry to be the same as overlayed widget
              this->setOverlayGeometry(widget->geometry(), false);
            }

            // Show and Register overlay
            this->Implementation->CheckOverlay->show();
            this->Implementation->CheckOverlayWidgetOn = widget;
          }
        }
        // Resize event
        if (event->type() == QEvent::Resize)
        {
          // Set overlay geometry
          QRect geometry = widget->geometry();
          if (this->Implementation->CheckOverlay->GlWidget)
          {
            geometry.moveTo(widget->mapToGlobal(QPoint(0, 0)));
          }
          this->setOverlayGeometry(geometry, false);
        }

        // Mouse button release -> Check Event
        if (event->type() == QEvent::MouseButtonRelease)
        {
          // Check Translators
          for (int i = 0; i != this->Implementation->Translators.size(); ++i)
          {
            bool error = false;
            if (this->Implementation->Translators[i]->translateEvent(
                  widget, event, pqEventTypes::CHECK_EVENT, error))
            {
              if (error)
              {
                qWarning() << "Error translating a check event for object " << widget;
              }
              return true;
            }
          }

          // Cannot check widget , Inform user trying to check uncheckable widget
          qWarning() << "Error checking an event for object, widget type not supported.";
          qWarning() << "Name of the widget:" << widget->objectName()
                     << ". Type of the widget:" << widget->metaObject()->className();
          if (widget->parent() != NULL)
          {
            qWarning() << "Type of parent widget:" << widget->parent()->metaObject()->className();
          }
        }
        // Block all input events, so the UI is static but still drawn.
        // Except for MouseMove
        if (dynamic_cast<QInputEvent*>(event) != NULL && event->type() != QEvent::MouseMove)
        {
          return true;
        }
      }
      else
      {
        // Event Recording
        for (int i = 0; i != this->Implementation->Translators.size(); ++i)
        {
          bool error = false;
          if (this->Implementation->Translators[i]->translateEvent(
                object, event, pqEventTypes::ACTION_EVENT, error))
          {
            if (error)
            {
              qWarning() << "Error translating an event for object " << object;
            }
            return false;
          }
        }
      }
    }
  }
  return false;
}

// ----------------------------------------------------------------------------
void pqEventTranslator::onRecordEvent(
  QObject* Object, const QString& Command, const QString& Arguments)
{
  this->onRecordEvent(Object, Command, Arguments, pqEventTypes::ACTION_EVENT);
}

// ----------------------------------------------------------------------------
void pqEventTranslator::onRecordEvent(
  QObject* Object, const QString& Command, const QString& Arguments, int eventType)
{
  if (this->Implementation->IgnoredObjects.contains(Object))
  {
    QRegularExpression commandFilter = this->Implementation->IgnoredObjects.value(Object);
    if (Command.contains(commandFilter))
    {
      return;
    }
  }

  QVariant blockRecordCommands = Object->property("BlockRecordCommands");
  if (blockRecordCommands.isValid() && Command.contains(blockRecordCommands.toRegularExpression()))
  {
    return;
  }

  QString name;
  if (eventType == pqEventTypes::ACTION_EVENT)
  {
    // When sender is pqEventObject, the Object name can be NULL.
    if (!qobject_cast<pqEventComment*>(this->sender()) || Object)
    {
      name = pqObjectNaming::GetName(*Object);
      if (name.isEmpty())
        return;
    }
  }
  else
  {
    // Check the QObject does have a name
    name = pqObjectNaming::GetName(*Object);
    if (name.isEmpty())
    {
      return;
    }
  }

  // Record user interaction time
  if (this->Implementation->RecordInteractionTimings)
  {
    if (this->Implementation->InteractionsTimer.isValid())
    {
      Q_EMIT recordEvent(name, "pause",
        QString::number(this->Implementation->InteractionsTimer.restart()),
        pqEventTypes::ACTION_EVENT);
    }
    else
    {
      this->Implementation->InteractionsTimer.start();
    }
  }

  // Record the event
  Q_EMIT recordEvent(name, Command, Arguments, eventType);
}

// ----------------------------------------------------------------------------
void pqEventTranslator::check(bool value)
{
  this->Implementation->Checking = value;

  // Hide overlay when not checking
  if (!value)
  {
    this->Implementation->hideOverlay();
  }
}

// ----------------------------------------------------------------------------
void pqEventTranslator::record(bool value)
{
  this->Implementation->Recording = value;

  // Hide overlay when not recording
  if (!value)
  {
    this->Implementation->hideOverlay();
  }
}

// ----------------------------------------------------------------------------
bool pqEventTranslator::isRecording()
{
  return this->Implementation->Recording;
}

// ----------------------------------------------------------------------------
void pqEventTranslator::recordInteractionTimings(bool value)
{
  if (value != this->Implementation->RecordInteractionTimings)
  {
    this->Implementation->RecordInteractionTimings = value;
    this->Implementation->InteractionsTimer.invalidate();
  }
}

// ----------------------------------------------------------------------------
void pqEventTranslator::setOverlayGeometry(const QRect& geometry, bool specific)
{
  if (this->Implementation->CheckOverlay != NULL)
  {
    this->Implementation->CheckOverlay->setGeometry(geometry);
  }
  this->Implementation->CheckOverlay->Specific = specific;
}

// ----------------------------------------------------------------------------
void pqEventTranslator::recordDashboardModeToggle(QObject* object, bool toggle)
{
  QString name = pqObjectNaming::GetName(*object);
  if (name.isEmpty())
  {
    qWarning() << "Error recording a dashboard mode event";
    return;
  }

  Q_EMIT recordEvent(
    name, "dashboard_mode", QVariant(toggle).toString(), pqEventTypes::ACTION_EVENT);
}
