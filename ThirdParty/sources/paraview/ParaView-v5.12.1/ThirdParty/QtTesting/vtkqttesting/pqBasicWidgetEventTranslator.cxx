/*=========================================================================

   Program: ParaView
   Module:    pqBasicWidgetEventTranslator.cxx

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

#include "pqBasicWidgetEventTranslator.h"

#include <QDialog>
#include <QEvent>
#include <QKeyEvent>
#include <QMetaProperty>
#include <QScrollBar>
#include <QWidget>

#include "pqEventTypes.h"

pqBasicWidgetEventTranslator::pqBasicWidgetEventTranslator(QObject* p)
  : pqWidgetEventTranslator(p)
{
}

pqBasicWidgetEventTranslator::~pqBasicWidgetEventTranslator()
{
}

bool pqBasicWidgetEventTranslator::translateEvent(
  QObject* object, QEvent* event, int eventType, bool& error)
{
  QWidget* widget = qobject_cast<QWidget*>(object);
  if (!widget)
    return false;

  if (eventType == pqEventTypes::ACTION_EVENT)
  {
    switch (event->type())
    {
      case QEvent::KeyPress:
      {
        QKeyEvent* keyEvent = static_cast<QKeyEvent*>(event);
        if (qobject_cast<QDialog*>(object))
        {
          Q_EMIT recordEvent(widget, "key", QString::number(keyEvent->key()));
        }
        return true;
        break;
      }
      case QEvent::MouseButtonPress:
      case QEvent::MouseButtonDblClick:
      case QEvent::MouseButtonRelease:
      {
        QMouseEvent* mouseEvent = static_cast<QMouseEvent*>(event);
        QString info = QString("%1,%2,%3,%4,%5")
                         .arg(mouseEvent->button())
                         .arg(mouseEvent->buttons())
                         .arg(mouseEvent->modifiers())
                         .arg(mouseEvent->x())
                         .arg(mouseEvent->y());

        if (event->type() != QEvent::MouseButtonRelease)
        {
          this->LastPos = mouseEvent->pos();
        }

        if (event->type() == QEvent::MouseButtonPress)
        {
          Q_EMIT recordEvent(widget, "mousePress", info);
        }
        if (event->type() == QEvent::MouseButtonDblClick)
        {
          Q_EMIT recordEvent(widget, "mouseDblClick", info);
        }
        else if (event->type() == QEvent::MouseButtonRelease)
        {
          if (this->LastPos != mouseEvent->pos())
          {
            Q_EMIT recordEvent(widget, "mouseMove", info);
          }
          Q_EMIT recordEvent(widget, "mouseRelease", info);
        }
        return true;
        break;
      }
      case QEvent::Wheel:
      {
        if (qobject_cast<QScrollBar*>(object))
        {
          QWheelEvent* wheelEvent = dynamic_cast<QWheelEvent*>(event);
          if (wheelEvent)
          {
            int buttons = wheelEvent->buttons();
            int modifiers = wheelEvent->modifiers();
            int numStep = wheelEvent->angleDelta().y();
            Q_EMIT recordEvent(object, "mouseWheel", QString("%1,%2,%3,%4,%5")
                                                       .arg(numStep)
                                                       .arg(buttons)
                                                       .arg(modifiers)
#if QT_VERSION >= QT_VERSION_CHECK(5, 15, 0)
                                                       .arg(wheelEvent->position().x())
                                                       .arg(wheelEvent->position().y()));
#else
                                                       .arg(wheelEvent->x())
                                                       .arg(wheelEvent->y()));
#endif
          }
        }
        return true;
        break;
      }
      default:
        break;
    }
  }
  else if (eventType == pqEventTypes::CHECK_EVENT)
  {
    if (event->type() == QEvent::MouseMove)
    {
      // Check for available meta prop
      const QMetaProperty metaProp = widget->metaObject()->userProperty();
      if (!metaProp.isValid() && qobject_cast<QWidget*>(widget->parent()) != NULL)
      {
        // MouseEvent can be received by the viewport
        const QMetaProperty metaProp = widget->parent()->metaObject()->userProperty();
      }
      if (metaProp.isValid())
      {
        return true;
      }
    }

    // Clicking while checking, actual check event
    if (event->type() == QEvent::MouseButtonRelease)
    {
      // Generic Meta prop check
      const QMetaProperty metaProp = widget->metaObject()->userProperty();
      if (!metaProp.isValid() && widget->parent() != NULL)
      {
        // MouseEvent can be received by the viewport, so try the parent widget
        const QMetaProperty metaProp = widget->parent()->metaObject()->userProperty();
        widget = qobject_cast<QWidget*>(widget->parent());
      }

      if (metaProp.isValid() && widget)
      {
        // Recover meto prop name
        QString propName = metaProp.name();

        // Record check event
        Q_EMIT recordEvent(widget, propName,
          widget->property(propName.toUtf8().data()).toString().replace("\t", " "),
          pqEventTypes::CHECK_EVENT);
        return true;
      }
    }
  }
  return this->Superclass::translateEvent(object, event, eventType, error);
}
