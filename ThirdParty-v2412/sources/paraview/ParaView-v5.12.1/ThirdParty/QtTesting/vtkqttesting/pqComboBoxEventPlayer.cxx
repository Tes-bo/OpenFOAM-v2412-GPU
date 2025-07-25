/*=========================================================================

  Program: ParaView
  Module:    pqComboBoxEventPlayer.cxx

  Copyright (c) 2005,2006 Sandia Corporation, Kitware Inc.
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

========================================================================*/
#include "pqComboBoxEventPlayer.h"
#include "pqEventTypes.h"
#include "pqObjectNaming.h"

#include <QComboBox>
#include <QDebug>

//-----------------------------------------------------------------------------
pqComboBoxEventPlayer::pqComboBoxEventPlayer(QObject* parentObject)
  : Superclass(parentObject)
{
}

//-----------------------------------------------------------------------------
pqComboBoxEventPlayer::~pqComboBoxEventPlayer()
{
}

//-----------------------------------------------------------------------------
bool pqComboBoxEventPlayer::playEvent(
  QObject* object, const QString& command, const QString& arguments, int eventType, bool& error)
{
  QComboBox* comboBox = qobject_cast<QComboBox*>(object);
  if (!comboBox)
  {
    return false;
  }

  if (eventType == pqEventTypes::ACTION_EVENT)
  {
    if (command == "activated" || command == "editTextChanged")
    {
      int index = comboBox->findText(arguments);
      if (index != -1)
      {
        comboBox->setCurrentIndex(index);
        if (command == "activated")
        {
          Q_EMIT comboBox->activated(index);
        }
      }
      else if (comboBox->isEditable() && command == "editTextChanged")
      {
        comboBox->setCurrentText(arguments);
      }
      else
      {
        QString possibles;
        for (int i = 0; i < comboBox->count(); i++)
        {
          possibles += QString("\t") + comboBox->itemText(i) + QString("\n");
        }
        qCritical() << "Unable to find " << arguments
                    << " in combo box: " << pqObjectNaming::GetName(*comboBox)
                    << "\nPossible values are:\n"
                    << possibles;
        error = true;
      }
    }
    return true;
  }
  return this->Superclass::playEvent(object, command, arguments, eventType, error);
}
