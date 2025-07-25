/*=========================================================================

   Program: ParaView
   Module:    pqTableViewEventTranslator.cxx

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
#include "pqTableViewEventTranslator.h"

#include <QHeaderView>
#include <QTableView>
//-----------------------------------------------------------------------------
pqTableViewEventTranslator::pqTableViewEventTranslator(QObject* parentObject)
  : Superclass(parentObject)
{
}

//-----------------------------------------------------------------------------
pqTableViewEventTranslator::~pqTableViewEventTranslator()
{
}

//-----------------------------------------------------------------------------
void pqTableViewEventTranslator::onEnteredCheck(const QModelIndex& item)
{
  // Recover Table View
  QTableView* tableView = qobject_cast<QTableView*>(this->AbstractItemView);

  // Get a visual rectangle of the item
  QRect visualRect = tableView->visualRect(item);

  // Translate the rect of margins and headers
  visualRect.translate(tableView->contentsMargins().left(), tableView->contentsMargins().top());
  visualRect.translate(
    tableView->verticalHeader()->width(), tableView->horizontalHeader()->height());

  // Stor item and signal that a specific overlay is ready to be drawn
  this->ModelItemCheck = &item;
  Q_EMIT this->specificOverlay(visualRect);
}

//-----------------------------------------------------------------------------
QAbstractItemView* pqTableViewEventTranslator::findCorrectedAbstractItemView(QObject* object) const
{
  // Ignore QHeaderView event specifically
  if (qobject_cast<QHeaderView*>(object))
  {
    return NULL;
  }

  QAbstractItemView* abstractItemView = qobject_cast<QTableView*>(object);
  if (!abstractItemView)
  {
    // mouse events go to the viewport widget
    abstractItemView = qobject_cast<QTableView*>(object->parent());
  }
  return abstractItemView;
}
