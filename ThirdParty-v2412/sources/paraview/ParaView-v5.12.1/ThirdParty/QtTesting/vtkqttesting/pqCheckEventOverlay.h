/*=========================================================================

   Program: ParaView
   Module:    pqEventTranslator.h

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

#ifndef _pqCheckEventOverlay_h
#define _pqCheckEventOverlay_h

#include <QWidget>
/**
Display a green or red overlay rectangle on top of parent widget
\sa pqEventTranslator
*/
class pqCheckEventOverlay : public QWidget
{
  Q_OBJECT
  typedef QWidget Superclass;

public:
  pqCheckEventOverlay(QWidget* parent = 0);

  // true if the overlayed widget can be checked, false otherwise
  bool Valid;

  // true if the overlayed widget in an opengl widget, false otherwise
  bool GlWidget;

  // true if the size of the overlay is specifically defined
  // via resize slot, false otherwise
  bool Specific;

  // Static overlauy margin
  static const int OVERLAY_MARGIN = 2;

  // Static overlay pen width
  static const int OVERLAY_PEN_WIDTH = 5;

protected:
  void paintEvent(QPaintEvent*) override;
};

#endif
