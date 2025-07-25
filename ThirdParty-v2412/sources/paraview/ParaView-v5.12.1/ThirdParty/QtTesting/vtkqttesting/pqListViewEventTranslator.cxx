/*=========================================================================

   Program: ParaView
   Module:    pqListViewEventTranslator.cxx

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
#include "pqListViewEventTranslator.h"

#include <QHeaderView>
#include <QListView>
//-----------------------------------------------------------------------------
pqListViewEventTranslator::pqListViewEventTranslator(QObject* parentObject)
  : Superclass(parentObject)
{
}

//-----------------------------------------------------------------------------
pqListViewEventTranslator::~pqListViewEventTranslator()
{
}

//-----------------------------------------------------------------------------
void pqListViewEventTranslator::onEnteredCheck(const QModelIndex& item)
{
  // Recover list View
  QListView* listView = qobject_cast<QListView*>(this->AbstractItemView);

  // Get a visual rectangle of the item
  QRect visualRect = listView->visualRect(item);

  // Translate the rect of margins and headers
  visualRect.translate(listView->contentsMargins().left(), listView->contentsMargins().top());

  // Stor item and signal that a specific overlay is ready to be drawn
  this->ModelItemCheck = &item;
  Q_EMIT this->specificOverlay(visualRect);
}

//-----------------------------------------------------------------------------
QAbstractItemView* pqListViewEventTranslator::findCorrectedAbstractItemView(QObject* object) const
{
  // Ignore QHeaderView event specifically
  if (qobject_cast<QHeaderView*>(object))
  {
    return NULL;
  }

  QAbstractItemView* abstractItemView = qobject_cast<QListView*>(object);
  if (!abstractItemView)
  {
    // mouse events go to the viewport widget
    abstractItemView = qobject_cast<QListView*>(object->parent());
  }
  return abstractItemView;
}
