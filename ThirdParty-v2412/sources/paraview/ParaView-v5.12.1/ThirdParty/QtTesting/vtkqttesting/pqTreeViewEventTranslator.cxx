/*=========================================================================

   Program: ParaView
   Module:    pqTreeViewEventTranslator.cxx

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

#include "pqTreeViewEventTranslator.h"

#include <QHeaderView>
#include <QTreeView>

//-----------------------------------------------------------------------------
pqTreeViewEventTranslator::pqTreeViewEventTranslator(QObject* parentObject)
  : Superclass(parentObject)
{
}

//-----------------------------------------------------------------------------
pqTreeViewEventTranslator::~pqTreeViewEventTranslator()
{
}

//-----------------------------------------------------------------------------
void pqTreeViewEventTranslator::monitorSignalsInternal(QAbstractItemView* abstractItemView)
{
  QObject::connect(abstractItemView, SIGNAL(expanded(const QModelIndex&)), this,
    SLOT(onExpanded(const QModelIndex&)));
  QObject::connect(abstractItemView, SIGNAL(collapsed(const QModelIndex&)), this,
    SLOT(onCollapsed(const QModelIndex&)));
}

//-----------------------------------------------------------------------------
void pqTreeViewEventTranslator::onExpanded(const QModelIndex& index)
{
  QTreeView* treeView = qobject_cast<QTreeView*>(this->sender());

  // record the check state change if the item is user-checkable.
  Q_EMIT this->recordEvent(treeView, "expand", this->getIndexAsString(index));
}

//-----------------------------------------------------------------------------
void pqTreeViewEventTranslator::onCollapsed(const QModelIndex& index)
{
  QTreeView* treeView = qobject_cast<QTreeView*>(this->sender());

  // record the check state change if the item is user-checkable.
  Q_EMIT this->recordEvent(treeView, "collapse", this->getIndexAsString(index));
}

//-----------------------------------------------------------------------------
void pqTreeViewEventTranslator::onEnteredCheck(const QModelIndex& item)
{
  // Recover treeview
  QTreeView* treeView = qobject_cast<QTreeView*>(this->AbstractItemView);

  // Recover a visual rectangle corresponding to the item
  QRect visualRect = treeView->visualRect(item);

  // Translating it of margins, as they are missing from visual rect
  visualRect.translate(treeView->contentsMargins().left(), treeView->contentsMargins().top());

  // As well as header
  visualRect.translate(0, treeView->header()->height());

  // store item and signal that a specific overlay is ready to be drawn
  this->ModelItemCheck = &item;
  Q_EMIT this->specificOverlay(visualRect);
}

//-----------------------------------------------------------------------------
QAbstractItemView* pqTreeViewEventTranslator::findCorrectedAbstractItemView(QObject* object) const
{
  // Ignore QHeaderView event specifically
  if (qobject_cast<QHeaderView*>(object))
  {
    return NULL;
  }

  QAbstractItemView* abstractItemView = qobject_cast<QTreeView*>(object);
  if (!abstractItemView)
  {
    // mouse events go to the viewport widget
    abstractItemView = qobject_cast<QTreeView*>(object->parent());
  }
  return abstractItemView;
}
