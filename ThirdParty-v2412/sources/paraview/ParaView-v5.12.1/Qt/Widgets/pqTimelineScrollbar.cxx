// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-FileCopyrightText: Copyright (c) Sandia Corporation
// SPDX-License-Identifier: BSD-3-Clause

// Hide PARAVIEW_DEPRECATED_IN_5_12_0() warnings for this class.
#define PARAVIEW_DEPRECATION_LEVEL 0

#include "pqTimelineScrollbar.h"

#include <QHBoxLayout>
#include <QScrollBar>

#include "pqAnimationModel.h"

constexpr int TimeScrollbarGranularity = 100000;

//-----------------------------------------------------------------------------
pqTimelineScrollbar::pqTimelineScrollbar(QWidget* p)
  : QWidget(p)
{
  QHBoxLayout* timeLayout = new QHBoxLayout(this);
  timeLayout->setContentsMargins(0, 0, 0, 0);
  this->TimeScrollBar = new QScrollBar;
  this->TimeScrollBar->setOrientation(Qt::Orientation::Horizontal);
  this->ScrollBarSpacer = new QSpacerItem(0, 0);
  timeLayout->addSpacerItem(this->ScrollBarSpacer);
  timeLayout->addWidget(this->TimeScrollBar);
  QObject::connect(this->TimeScrollBar, SIGNAL(actionTriggered(int)), this, SLOT(setTimeZoom(int)));
}

//-----------------------------------------------------------------------------
void pqTimelineScrollbar::setAnimationModel(pqAnimationModel* model)
{
  if (this->AnimationModel)
  {
    this->AnimationModel->disconnect(this);
  }

  this->AnimationModel = model;

  if (!this->AnimationModel)
  {
    return;
  }

  QObject::connect(this->AnimationModel, SIGNAL(zoomChanged()), this, SLOT(updateTimeScrollbar()));
  this->updateTimeScrollbar();
}

//-----------------------------------------------------------------------------
void pqTimelineScrollbar::linkSpacing(QObject* spaceNotifier)
{

  if (this->SpacingNotifier)
  {
    this->AnimationModel->disconnect(this);
  }

  this->SpacingNotifier = spaceNotifier;

  if (!this->SpacingNotifier)
  {
    return;
  }

  QObject::connect(
    spaceNotifier, SIGNAL(timelineOffsetChanged(int)), this, SLOT(updateTimeScrollbarOffset(int)));
}

//-----------------------------------------------------------------------------
void pqTimelineScrollbar::updateTimeScrollbarOffset(int offset)
{
  this->ScrollBarSpacer->changeSize(offset, 0);
}

//-----------------------------------------------------------------------------
void pqTimelineScrollbar::setTimeZoom(int action)
{
  pqAnimationModel* animModel = this->AnimationModel;

  if (!animModel)
  {
    return;
  }

  double timeValue;

  switch (action)
  {
    case QAbstractSlider::SliderSingleStepAdd:
      timeValue = this->TimeScrollBar->value() + this->TimeScrollBar->pageStep();
      this->TimeScrollBar->setValue(timeValue);
      break;
    case QAbstractSlider::SliderSingleStepSub:
      timeValue = this->TimeScrollBar->value() - this->TimeScrollBar->pageStep();
      this->TimeScrollBar->setValue(timeValue);
      break;
    default:
      timeValue = this->TimeScrollBar->sliderPosition();
  }

  if (timeValue >= this->TimeScrollBar->maximum())
  {
    timeValue = animModel->endTime() -
      (animModel->endTime() - animModel->startTime()) / animModel->zoomFactor();
  }
  else if (timeValue <= this->TimeScrollBar->minimum())
  {
    timeValue = animModel->startTime();
  }
  else
  {
    timeValue =
      (timeValue / TimeScrollbarGranularity) * (animModel->endTime() - animModel->startTime()) +
      animModel->startTime();
  }
  animModel->positionZoom(timeValue);
}

//-----------------------------------------------------------------------------
void pqTimelineScrollbar::updateTimeScrollbar()
{
  pqAnimationModel* animModel = this->AnimationModel;

  if (!animModel)
  {
    return;
  }

  QScrollBar* scrollb = this->TimeScrollBar;

  if (animModel->zoomFactor() == 1)
  {
    scrollb->setVisible(false);
    return;
  }

  scrollb->setVisible(true);
  int rangeMax = TimeScrollbarGranularity - (TimeScrollbarGranularity / animModel->zoomFactor());
  scrollb->setPageStep(TimeScrollbarGranularity / animModel->zoomFactor());
  scrollb->setRange(0, rangeMax);
  int relativeTime = ((animModel->zoomStartTime() - animModel->startTime()) /
                       (animModel->endTime() - animModel->startTime())) *
    TimeScrollbarGranularity;
  scrollb->setValue(relativeTime);
}
