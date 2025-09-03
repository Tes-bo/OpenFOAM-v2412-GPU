/*=========================================================================

   Program: ParaView
   Module:    pqLineEditEventTranslator.cxx

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

#include "pqLineEditEventTranslator.h"

#include <QDebug>
#include <QEvent>
#include <QKeyEvent>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QSpinBox>
#include <QTextDocument>
#include <QTextEdit>

#include "pqEventTypes.h"

pqLineEditEventTranslator::pqLineEditEventTranslator(QObject* p)
  : pqWidgetEventTranslator(p)
{
}

bool pqLineEditEventTranslator::translateEvent(
  QObject* object, QEvent* event, int eventType, bool& error)
{
  QObject* tmpObject = object;
  QLineEdit* leObject = qobject_cast<QLineEdit*>(object);
  QTextEdit* teObject = qobject_cast<QTextEdit*>(object);
  QPlainTextEdit* pteObject = qobject_cast<QPlainTextEdit*>(object);
  if (!leObject && !teObject && !pteObject && object->parent() != NULL)
  {
    // MouseEvent can be received by viewport
    tmpObject = object->parent();
    leObject = qobject_cast<QLineEdit*>(tmpObject);
    teObject = qobject_cast<QTextEdit*>(tmpObject);
    pteObject = qobject_cast<QPlainTextEdit*>(tmpObject);
  }
  if (!leObject && !teObject && !pteObject)
  {
    return false;
  }

  if (eventType == pqEventTypes::ACTION_EVENT)
  {
    // If this line edit is part of a spinbox, don't translate events
    // (the spinbox translator will receive the final value directly)
    if (qobject_cast<QSpinBox*>(tmpObject->parent()))
    {
      return false;
    }

    switch (event->type())
    {
      case QEvent::KeyRelease:
      {
        QKeyEvent* ke = static_cast<QKeyEvent*>(event);
        QString keyText = ke->text();
        if (keyText.length() && keyText.at(0).isPrint())
        {
          if (leObject)
          {
            Q_EMIT recordEvent(tmpObject, "set_string", leObject->text());
          }
          else if (teObject)
          {
            Q_EMIT recordEvent(tmpObject, "set_string", teObject->document()->toPlainText());
          }
          else if (pteObject)
          {
            Q_EMIT recordEvent(tmpObject, "set_string", pteObject->document()->toPlainText());
          }
        }
        // if we record F2 event, will cause some issue with the TreeView
        // Need test to know if we need to record those events
        else if (ke->key() != Qt::Key_F2)
        {
          Q_EMIT recordEvent(tmpObject, "key", QString("%1").arg(ke->key()));
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
    // TextEdit only, LineEdit does not need specific check
    if (teObject != NULL || pteObject != NULL)
    {
      if (event->type() == QEvent::MouseMove)
      {
        return true;
      }
      // Clicking while checking, actual check event
      if (event->type() == QEvent::MouseButtonRelease)
      {
        if (teObject != NULL)
        {
          Q_EMIT this->recordEvent(teObject, "plainText",
            teObject->toPlainText().replace("\t", " "), pqEventTypes::CHECK_EVENT);
        }
        else /* if (pteObject != NULL)*/
        {
          Q_EMIT this->recordEvent(pteObject, "plainText",
            pteObject->toPlainText().replace("\t", " "), pqEventTypes::CHECK_EVENT);
        }
        return true;
      }
    }
  }
  return this->Superclass::translateEvent(object, event, eventType, error);
}
