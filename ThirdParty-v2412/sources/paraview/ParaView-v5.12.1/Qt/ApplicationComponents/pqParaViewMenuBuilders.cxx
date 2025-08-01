// SPDX-FileCopyrightText: Copyright (c) Kitware Inc.
// SPDX-FileCopyrightText: Copyright (c) Sandia Corporation
// SPDX-License-Identifier: BSD-3-Clause
#include "pqParaViewMenuBuilders.h"

#include "ui_pqEditMenuBuilder.h"
#include "ui_pqFileMenuBuilder.h"

#include "pqQtConfig.h"

#include "pqAboutDialogReaction.h"
#include "pqAnimatedExportReaction.h"
#include "pqAnimationTimeToolbar.h"
#include "pqApplicationCore.h"
#include "pqApplicationSettingsReaction.h"
#include "pqApplyPropertiesReaction.h"
#include "pqAxesToolbar.h"
#include "pqCameraLinkReaction.h"
#include "pqCameraToolbar.h"
#include "pqCameraUndoRedoReaction.h"
#include "pqCatalystConnectReaction.h"
#include "pqCatalystContinueReaction.h"
#include "pqCatalystExportReaction.h"
#include "pqCatalystPauseSimulationReaction.h"
#include "pqCatalystRemoveBreakpointReaction.h"
#include "pqCatalystSetBreakpointReaction.h"
#include "pqCategoryToolbarsBehavior.h"
#include "pqChangeFileNameReaction.h"
#include "pqChangePipelineInputReaction.h"
#include "pqColorToolbar.h"
#include "pqCopyReaction.h"
#include "pqCreateCustomFilterReaction.h"
#include "pqCustomViewpointsDefaultController.h"
#include "pqCustomViewpointsToolbar.h"
#include "pqCustomizeShortcutsReaction.h"
#include "pqDataQueryReaction.h"
#include "pqDeleteReaction.h"
#include "pqDesktopServicesReaction.h"
#include "pqExampleVisualizationsDialogReaction.h"
#include "pqExportReaction.h"
#include "pqExtractorsMenuReaction.h"
#include "pqFiltersMenuReaction.h"
#ifdef PARAVIEW_USE_QTHELP
#include "pqHelpReaction.h"
#endif
#include "pqIgnoreSourceTimeReaction.h"
#include "pqLightToolbar.h"
#include "pqLinkSelectionReaction.h"
#include "pqLoadDataReaction.h"
#include "pqLoadMaterialsReaction.h"
#include "pqLoadRestoreWindowLayoutReaction.h"
#include "pqLoadStateReaction.h"
#include "pqLogViewerReaction.h"
#include "pqMainControlsToolbar.h"
#include "pqManageCustomFiltersReaction.h"
#include "pqManageExpressionsReaction.h"
#include "pqManageFavoritesReaction.h"
#include "pqManageLinksReaction.h"
#include "pqManagePluginsReaction.h"
#include "pqPVApplicationCore.h"
#include "pqPropertiesPanel.h"
#include "pqProxyGroupMenuManager.h"
#include "pqRecentFilesMenu.h"
#include "pqReloadFilesReaction.h"
#include "pqRenameProxyReaction.h"
#include "pqRepresentationToolbar.h"
#include "pqResetDefaultSettingsReaction.h"
#include "pqSaveAnimationGeometryReaction.h"
#include "pqSaveAnimationReaction.h"
#include "pqSaveDataReaction.h"
#include "pqSaveExtractsReaction.h"
#include "pqSaveScreenshotReaction.h"
#include "pqSaveStateReaction.h"
#include "pqSearchItemReaction.h"
#include "pqServerConnectReaction.h"
#include "pqServerDisconnectReaction.h"
#include "pqSetMainWindowTitleReaction.h"
#include "pqSetName.h"
#include "pqShowHideAllReaction.h"
#include "pqSourcesMenuReaction.h"
#include "pqTemporalExportReaction.h"
#include "pqTestingReaction.h"
#include "pqTimerLogReaction.h"
#include "pqUndoRedoReaction.h"
#include "pqVCRToolbar.h"
#include "pqViewMenuManager.h"

#if VTK_MODULE_ENABLE_ParaView_pqPython
#include "pqEditMacrosReaction.h"
#include "pqMacroReaction.h"
#include "pqPythonMacroSupervisor.h"
#include "pqPythonManager.h"
#include "pqPythonScriptEditorReaction.h"
#include "pqPythonTabWidget.h"
#include "pqTraceReaction.h"
#endif

#include <QApplication>
#include <QCoreApplication>
#include <QDockWidget>
#include <QFileInfo>
#include <QKeySequence>
#include <QLayout>
#include <QMainWindow>
#include <QMenu>
#include <QProxyStyle>
#include <QSysInfo>

#include "vtkPVFileInformation.h"
#include "vtkPVVersion.h"
#include "vtkSMProxyManager.h"

//-----------------------------------------------------------------------------
class pqActiveDisabledStyle : public QProxyStyle
{
public:
  int styleHint(StyleHint hint, const QStyleOption* option = nullptr,
    const QWidget* widget = nullptr, QStyleHintReturn* returnData = nullptr) const override
  {
    return hint == QStyle::SH_Menu_AllowActiveAndDisabled
      ? 1
      : QProxyStyle::styleHint(hint, option, widget, returnData);
  }
};

//-----------------------------------------------------------------------------
void pqParaViewMenuBuilders::buildFileMenu(QMenu& menu)
{
  QString objectName = menu.objectName();
  Ui::pqFileMenuBuilder ui;
  ui.setupUi(&menu);
  // since the UI file tends to change the name of the menu.
  menu.setObjectName(objectName);

  QObject::connect(
    ui.actionFileExit, SIGNAL(triggered()), QApplication::instance(), SLOT(closeAllWindows()));

  // now setup reactions.
  new pqLoadDataReaction(ui.actionFileOpen);
#if VTK_MODULE_ENABLE_VTK_RenderingRayTracing
  new pqLoadMaterialsReaction(ui.actionFileLoadMaterials);
#else
  ui.actionFileLoadMaterials->setEnabled(false);
#endif
  new pqRecentFilesMenu(*ui.menuRecentFiles, ui.menuRecentFiles);
  new pqReloadFilesReaction(ui.actionReloadFiles);

  new pqLoadStateReaction(ui.actionFileLoadServerState);
  new pqSaveStateReaction(ui.actionFileSaveServerState);
#if VTK_MODULE_ENABLE_ParaView_pqPython
  new pqCatalystExportReaction(ui.actionFileSaveCatalystState);
#else
  ui.actionFileSaveCatalystState->setEnabled(false);
#endif

  new pqServerConnectReaction(ui.actionServerConnect);
  new pqServerDisconnectReaction(ui.actionServerDisconnect);

  new pqSaveScreenshotReaction(ui.actionFileSaveScreenshot);
  new pqSaveAnimationReaction(ui.actionFileSaveAnimation);
  new pqSaveAnimationGeometryReaction(ui.actionFileSaveGeometry);

  new pqExportReaction(ui.actionExport);
#if VTK_MODULE_ENABLE_ParaView_pqPython
  new pqAnimatedExportReaction(ui.actionAnimatedExport);
#else
  ui.actionAnimatedExport->setEnabled(false);
#endif
  new pqSaveExtractsReaction(ui.actionFileSaveExtracts);
  new pqSaveDataReaction(ui.actionFileSaveData);

  new pqLoadRestoreWindowLayoutReaction(true, ui.actionFileLoadWindowArrangement);
  new pqLoadRestoreWindowLayoutReaction(false, ui.actionFileSaveWindowArrangement);
}

//-----------------------------------------------------------------------------
void pqParaViewMenuBuilders::buildEditMenu(QMenu& menu, pqPropertiesPanel* propertiesPanel)
{
  QString objectName = menu.objectName();
  Ui::pqEditMenuBuilder ui;
  ui.setupUi(&menu);
  // since the UI file tends to change the name of the menu.
  menu.setObjectName(objectName);

  new pqUndoRedoReaction(ui.actionEditUndo, true);
  new pqUndoRedoReaction(ui.actionEditRedo, false);
  new pqCameraUndoRedoReaction(ui.actionEditCameraUndo, true);
  new pqCameraUndoRedoReaction(ui.actionEditCameraRedo, false);
  new pqChangePipelineInputReaction(ui.actionChangeInput);
  new pqChangeFileNameReaction(ui.actionChangeFile);
  new pqIgnoreSourceTimeReaction(ui.actionIgnoreTime);
  new pqDeleteReaction(ui.actionDelete);
  ui.actionDelete->setShortcut(QKeySequence(Qt::ALT + Qt::Key_D));
  new pqDeleteReaction(ui.actionDeleteTree, pqDeleteReaction::DeleteModes::TREE);
  new pqDeleteReaction(ui.actionDelete_All, pqDeleteReaction::DeleteModes::ALL);
  new pqShowHideAllReaction(ui.actionShow_All, pqShowHideAllReaction::ActionType::Show);
  new pqShowHideAllReaction(ui.actionHide_All, pqShowHideAllReaction::ActionType::Hide);
  new pqSaveScreenshotReaction(ui.actionCopyScreenshotToClipboard, true);
  new pqCopyReaction(ui.actionCopy);
  new pqCopyReaction(ui.actionPaste, true);
  new pqCopyReaction(ui.actionCopyPipeline, false, true);
  new pqCopyReaction(ui.actionPastePipeline, true, true);
  new pqApplicationSettingsReaction(ui.actionEditSettings);
  new pqDataQueryReaction(ui.actionQuery);
  new pqSearchItemReaction(ui.actionSearch);
  new pqResetDefaultSettingsReaction(ui.actionResetDefaultSettings);
  new pqSetMainWindowTitleReaction(ui.actionSetMainWindowTitle);
  new pqRenameProxyReaction(ui.actionRename, propertiesPanel);

  if (propertiesPanel)
  {
    QAction* applyAction = new QAction(QIcon(":/pqWidgets/Icons/pqApply.svg"),
      QCoreApplication::translate("pqEditMenu", "Apply"), &menu);
    applyAction->setShortcut(QKeySequence(Qt::ALT + Qt::Key_A));
    QAction* resetAction = new QAction(QIcon(":/pqWidgets/Icons/pqCancel.svg"),
      QCoreApplication::translate("pqEditMenu", "Reset"), &menu);
    resetAction->setShortcut(QKeySequence(Qt::ALT + Qt::Key_R));
    menu.insertAction(ui.actionDelete, applyAction);
    menu.insertAction(ui.actionDelete, resetAction);
    new pqApplyPropertiesReaction(propertiesPanel, applyAction, true);
    new pqApplyPropertiesReaction(propertiesPanel, resetAction, false);
  }
}

//-----------------------------------------------------------------------------
void pqParaViewMenuBuilders::buildSourcesMenu(QMenu& menu, QMainWindow* mainWindow)
{
  pqProxyGroupMenuManager* mgr = new pqProxyGroupMenuManager(&menu, "ParaViewSources");
  mgr->addProxyDefinitionUpdateListener("sources");
  mgr->setRecentlyUsedMenuSize(10);
  new pqSourcesMenuReaction(mgr);
  if (mainWindow)
  {
    // create toolbars for categories as needed.
    new pqCategoryToolbarsBehavior(mgr, mainWindow);
  }
}

//-----------------------------------------------------------------------------
void pqParaViewMenuBuilders::buildFiltersMenu(
  QMenu& menu, QMainWindow* mainWindow, bool hideDisabled, bool quickLaunchable)
{
  // Make sure disabled actions are still considered active
  auto* style = new pqActiveDisabledStyle;
  style->setParent(&menu);
  menu.setStyle(style);

  pqProxyGroupMenuManager* mgr =
    new pqProxyGroupMenuManager(&menu, "ParaViewFilters", quickLaunchable);
  mgr->addProxyDefinitionUpdateListener("filters");
  mgr->setRecentlyUsedMenuSize(10);
  mgr->setEnableFavorites(true);
  pqFiltersMenuReaction* menuReaction = new pqFiltersMenuReaction(mgr, hideDisabled);

  // Connect the filters menu about to show and the quick-launch dialog about to show
  // to update the enabled/disabled state of the menu items via the
  // pqFiltersMenuReaction
  menuReaction->connect(&menu, SIGNAL(aboutToShow()), SLOT(updateEnableState()));

  auto* pvappcore = pqPVApplicationCore::instance();
  if (quickLaunchable && pvappcore)
  {
    menuReaction->connect(pvappcore, SIGNAL(aboutToShowQuickLaunch()), SLOT(updateEnableState()));
  }

  if (mainWindow)
  {
    // create toolbars for categories as needed.
    new pqCategoryToolbarsBehavior(mgr, mainWindow);
  }
}

//-----------------------------------------------------------------------------
void pqParaViewMenuBuilders::buildExtractorsMenu(
  QMenu& menu, QMainWindow* mainWindow, bool hideDisabled, bool quickLaunchable)
{
  // Make sure disabled actions are still considered active
  auto* style = new pqActiveDisabledStyle;
  style->setParent(&menu);
  menu.setStyle(style);

  pqProxyGroupMenuManager* mgr =
    new pqProxyGroupMenuManager(&menu, "ParaViewExtractWriters", quickLaunchable);
  mgr->addProxyDefinitionUpdateListener("extract_writers");

  auto menuReaction = new pqExtractorsMenuReaction(mgr, hideDisabled);
  auto* pvappcore = pqPVApplicationCore::instance();
  if (quickLaunchable && pvappcore)
  {
    menuReaction->connect(pvappcore, SIGNAL(aboutToShowQuickLaunch()), SLOT(updateEnableState()));
  }

  if (mainWindow)
  {
    // create toolbars for categories as needed.
    new pqCategoryToolbarsBehavior(mgr, mainWindow);
  }
}

//-----------------------------------------------------------------------------
void pqParaViewMenuBuilders::buildToolsMenu(QMenu& menu)
{
  new pqCreateCustomFilterReaction(
    menu.addAction(QCoreApplication::translate("pqToolsMenu", "Create Custom Filter..."))
    << pqSetName("actionToolsCreateCustomFilter"));
  new pqCameraLinkReaction(
    menu.addAction(QCoreApplication::translate("pqToolsMenu", "Add Camera Link..."))
    << pqSetName("actionToolsAddCameraLink"));
  new pqLinkSelectionReaction(
    menu.addAction(QCoreApplication::translate("pqToolsMenu", "Link with Selection"))
    << pqSetName("actionToolsLinkSelection"));
  menu.addSeparator();
  new pqManageCustomFiltersReaction(
    menu.addAction(QCoreApplication::translate("pqToolsMenu", "Manage Custom Filters..."))
    << pqSetName("actionToolsManageCustomFilters"));
  new pqManageLinksReaction(
    menu.addAction(QCoreApplication::translate("pqToolsMenu", "Manage Links..."))
    << pqSetName("actionToolsManageLinks"));
  //<addaction name="actionToolsAddCameraLink" />
  // Add support for importing plugins only if ParaView was built shared.
  new pqManagePluginsReaction(
    menu.addAction(QCoreApplication::translate("pqToolsMenu", "Manage Plugins..."))
    << pqSetName("actionManage_Plugins"));

  QMenu* dummyMenu = new QMenu();
  pqProxyGroupMenuManager* mgr = new pqProxyGroupMenuManager(dummyMenu, "ParaViewFilters", false);
  mgr->addProxyDefinitionUpdateListener("filters");

  QAction* manageFavoritesAction =
    menu.addAction(QCoreApplication::translate("pqToolsMenu", "Manage Favorites..."))
    << pqSetName("actionManage_Favorites");
  new pqManageFavoritesReaction(manageFavoritesAction, mgr);

  new pqCustomizeShortcutsReaction(
    menu.addAction(QCoreApplication::translate("pqToolsMenu", "Customize Shortcuts..."))
    << pqSetName("actionCustomize"));

  new pqManageExpressionsReaction(
    menu.addAction(QCoreApplication::translate("pqToolsMenu", "Manage Expressions"))
    << pqSetName("actionToolsManageExpressions"));

  menu.addSeparator(); // --------------------------------------------------

  //<addaction name="actionToolsDumpWidgetNames" />
  new pqTestingReaction(menu.addAction(QCoreApplication::translate("pqToolsMenu", "Record Test..."))
      << pqSetName("actionToolsRecordTest"),
    pqTestingReaction::RECORD);
  new pqTestingReaction(menu.addAction(QCoreApplication::translate("pqToolsMenu", "Play Test..."))
      << pqSetName("actionToolsPlayTest"),
    pqTestingReaction::PLAYBACK, Qt::QueuedConnection);
  new pqTestingReaction(menu.addAction(QCoreApplication::translate("pqToolsMenu", "Lock View Size"))
      << pqSetName("actionTesting_Window_Size"),
    pqTestingReaction::LOCK_VIEW_SIZE);
  new pqTestingReaction(
    menu.addAction(QCoreApplication::translate("pqToolsMenu", "Lock View Size Custom..."))
      << pqSetName("actionTesting_Window_Size_Custom"),
    pqTestingReaction::LOCK_VIEW_SIZE_CUSTOM);
  menu.addSeparator();
  new pqTimerLogReaction(menu.addAction(QCoreApplication::translate("pqToolsMenu", "Timer Log"))
    << pqSetName("actionToolsTimerLog"));
  new pqLogViewerReaction(menu.addAction(QCoreApplication::translate("pqToolsMenu", "Log Viewer"))
    << pqSetName("actionToolsLogViewer"));

#if VTK_MODULE_ENABLE_ParaView_pqPython
  menu.addSeparator(); // --------------------------------------------------
  new pqTraceReaction(menu.addAction(QCoreApplication::translate("pqToolsMenu", "Start Trace"))
      << pqSetName("actionToolsStartStopTrace"),
    QCoreApplication::translate("pqToolsMenu", "Start Trace"),
    QCoreApplication::translate("pqToolsMenu", "Stop Trace"));
  menu.addSeparator();
  new pqPythonScriptEditorReaction(
    menu.addAction(QCoreApplication::translate("pqToolsMenu", "Python Script Editor"))
    << pqSetName("actionToolsOpenScriptEditor"));
#endif
}

//-----------------------------------------------------------------------------
void pqParaViewMenuBuilders::buildViewMenu(QMenu& menu, QMainWindow& mainWindow)
{
  new pqViewMenuManager(&mainWindow, &menu);
}

//-----------------------------------------------------------------------------
void pqParaViewMenuBuilders::buildPipelineBrowserContextMenu(QMenu& menu, QMainWindow* mainWindow)
{
  // Build the context menu manually so we can insert submenus where needed.
  QAction* actionPBOpen = new QAction(menu.parent());
  actionPBOpen->setObjectName(QStringLiteral("actionPBOpen"));
  QIcon icon4;
  icon4.addFile(QStringLiteral(":/pqWidgets/Icons/pqOpen.svg"), QSize(), QIcon::Normal, QIcon::Off);
  actionPBOpen->setIcon(icon4);
  actionPBOpen->setShortcutContext(Qt::WidgetShortcut);
  actionPBOpen->setText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "&Open", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
  actionPBOpen->setToolTip(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Open", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_STATUSTIP
  actionPBOpen->setStatusTip(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Open", Q_NULLPTR));
#endif // QT_NO_STATUSTIP

  QAction* actionPBShowAll = new QAction(menu.parent());
  actionPBShowAll->setObjectName(QStringLiteral("actionPBShowAll"));
  QIcon showAllIcon;
  showAllIcon.addFile(
    QStringLiteral(":/pqWidgets/Icons/pqEyeball.svg"), QSize(), QIcon::Normal, QIcon::Off);
  actionPBShowAll->setIcon(showAllIcon);
  actionPBShowAll->setText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "&Show All", Q_NULLPTR));
#ifndef QT_NO_STATUSTIP
  actionPBShowAll->setStatusTip(QCoreApplication::translate(
    "pqPipelineBrowserContextMenu", "Show all source outputs in the pipeline", Q_NULLPTR));
#endif // QT_NO_STATUSTIP

  QAction* actionPBHideAll = new QAction(menu.parent());
  actionPBHideAll->setObjectName(QStringLiteral("actionPBHideAll"));
  QIcon hideAllIcon;
  hideAllIcon.addFile(
    QStringLiteral(":/pqWidgets/Icons/pqEyeballClosed.svg"), QSize(), QIcon::Normal, QIcon::Off);
  actionPBHideAll->setIcon(hideAllIcon);
  actionPBHideAll->setText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "&Hide All", Q_NULLPTR));
#ifndef QT_NO_STATUSTIP
  actionPBHideAll->setStatusTip(QCoreApplication::translate(
    "pqPipelineBrowserContextMenu", "Hide all source outputs in the pipeline", Q_NULLPTR));
#endif // QT_NO_STATUSTIP

  QAction* actionPBCopy = new QAction(menu.parent());
  actionPBCopy->setObjectName(QStringLiteral("actionPBCopy"));
  QIcon icon2;
  icon2.addFile(QStringLiteral(":/pqWidgets/Icons/pqCopy.svg"), QSize(), QIcon::Normal, QIcon::Off);
  actionPBCopy->setIcon(icon2);
  actionPBCopy->setText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "&Copy", Q_NULLPTR));
#ifndef QT_NO_STATUSTIP
  actionPBCopy->setStatusTip(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Copy Properties", Q_NULLPTR));
#endif // QT_NO_STATUSTIP

  QAction* actionPBPaste = new QAction(menu.parent());
  actionPBPaste->setObjectName(QStringLiteral("actionPBPaste"));
  QIcon icon3;
  icon3.addFile(
    QStringLiteral(":/pqWidgets/Icons/pqPaste.svg"), QSize(), QIcon::Normal, QIcon::Off);
  actionPBPaste->setIcon(icon3);
  actionPBPaste->setText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "&Paste", Q_NULLPTR));
#ifndef QT_NO_STATUSTIP
  actionPBPaste->setStatusTip(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Paste Properties", Q_NULLPTR));
#endif // QT_NO_STATUSTIP

  QAction* actionPBCopyPipeline = new QAction(menu.parent());
  actionPBCopyPipeline->setObjectName(QStringLiteral("actionPBCopyPipeline"));
  actionPBCopyPipeline->setIcon(icon2);
  actionPBCopyPipeline->setText(
    QApplication::translate("pqPipelineBrowserContextMenu", "Copy Pipeline", Q_NULLPTR));
#ifndef QT_NO_STATUSTIP
  actionPBCopyPipeline->setStatusTip(
    QApplication::translate("pqPipelineBrowserContextMenu", "Copy Pipeline", Q_NULLPTR));
#endif // QT_NO_STATUSTIP

  QAction* actionPBPastePipeline = new QAction(menu.parent());
  actionPBPastePipeline->setObjectName(QStringLiteral("actionPBPastePipeline"));
  actionPBPastePipeline->setIcon(icon3);
  actionPBPastePipeline->setText(
    QApplication::translate("pqPipelineBrowserContextMenu", "Paste Pipeline", Q_NULLPTR));
#ifndef QT_NO_STATUSTIP
  actionPBPastePipeline->setStatusTip(
    QApplication::translate("pqPipelineBrowserContextMenu", "Paste Pipeline", Q_NULLPTR));
#endif // QT_NO_STATUSTIP

  QAction* actionPBChangeInput = new QAction(menu.parent());
  actionPBChangeInput->setObjectName(QStringLiteral("actionPBChangeInput"));
  actionPBChangeInput->setText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Change &Input...", Q_NULLPTR));
  actionPBChangeInput->setIconText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Change Input...", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
  actionPBChangeInput->setToolTip(QCoreApplication::translate(
    "pqPipelineBrowserContextMenu", "Change a Filter's Input", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_STATUSTIP
  actionPBChangeInput->setStatusTip(QCoreApplication::translate(
    "pqPipelineBrowserContextMenu", "Change a Filter's Input", Q_NULLPTR));
#endif // QT_NO_STATUSTIP

  QAction* actionPBReloadFiles = new QAction(menu.parent());
  actionPBReloadFiles->setObjectName(QStringLiteral("actionPBReloadFiles"));
  actionPBReloadFiles->setText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Reload Files", Q_NULLPTR));
  actionPBReloadFiles->setIconText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Reload Files", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
  actionPBReloadFiles->setToolTip(QCoreApplication::translate("pqPipelineBrowserContextMenu",
    "Reload data files in case they were changed externally.", Q_NULLPTR));
#endif // QT_NO_TOOLTIP

  QAction* actionPBChangeFile = new QAction(menu.parent());
  actionPBChangeFile->setObjectName(QStringLiteral("actionPBChangeFile"));
  actionPBChangeFile->setText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Change File", Q_NULLPTR));
#ifndef QT_NO_STATUSTIP
  actionPBChangeFile->setStatusTip(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Change File", Q_NULLPTR));
#endif // QT_NO_STATUSTIP

  QAction* actionPBIgnoreTime = new QAction(menu.parent());
  actionPBIgnoreTime->setObjectName(QStringLiteral("actionPBIgnoreTime"));
  actionPBIgnoreTime->setCheckable(true);
  actionPBIgnoreTime->setText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Ignore Time", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
  actionPBIgnoreTime->setToolTip(QCoreApplication::translate("pqPipelineBrowserContextMenu",
    "Disregard this source/filter's time from animations", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_STATUSTIP
  actionPBIgnoreTime->setStatusTip(QCoreApplication::translate("pqPipelineBrowserContextMenu",
    "Disregard this source/filter's time from animations", Q_NULLPTR));
#endif // QT_NO_STATUSTIP

  QAction* actionPBDelete = new QAction(menu.parent());
  actionPBDelete->setObjectName(QStringLiteral("actionPBDelete"));
  QIcon icon;
  icon.addFile(
    QStringLiteral(":/QtWidgets/Icons/pqDelete.svg"), QSize(), QIcon::Normal, QIcon::Off);
  actionPBDelete->setIcon(icon);
  actionPBDelete->setText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "&Delete", Q_NULLPTR));
#ifndef QT_NO_STATUSTIP
  actionPBDelete->setStatusTip(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Delete", Q_NULLPTR));
#endif // QT_NO_STATUSTIP

  QByteArray signalName = QMetaObject::normalizedSignature("deleteKey()");
  if (menu.parent()->metaObject()->indexOfSignal(signalName) != -1)
  {
    // Trigger a delete when the user requests a delete.
    QObject::connect(
      menu.parent(), SIGNAL(deleteKey()), actionPBDelete, SLOT(trigger()), Qt::QueuedConnection);
  }

  QAction* actionPBDeleteTree = new QAction(menu.parent());
  actionPBDeleteTree->setObjectName(QStringLiteral("actionPBDeleteTree"));
  actionPBDeleteTree->setIcon(icon);
  actionPBDeleteTree->setText(QCoreApplication::translate(
    "pqPipelineBrowserContextMenu", "Delete Downstream Pipeline", Q_NULLPTR));
#ifndef QT_NO_STATUSTIP
  actionPBDeleteTree->setStatusTip(QCoreApplication::translate(
    "pqPipelineBrowserContextMenu", "Delete selection and all downstream filters", Q_NULLPTR));
#endif // QT_NO_STATUSTIP

  QAction* actionPBCreateCustomFilter = new QAction(menu.parent());
  actionPBCreateCustomFilter->setObjectName(QStringLiteral("actionPBCreateCustomFilter"));
  actionPBCreateCustomFilter->setText(QCoreApplication::translate(
    "pqPipelineBrowserContextMenu", "&Create Custom Filter...", Q_NULLPTR));

  QAction* actionPBLinkSelection = new QAction(menu.parent());
  actionPBLinkSelection->setObjectName(QStringLiteral("actionPBLinkSelection"));
  actionPBLinkSelection->setText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Link with selection", Q_NULLPTR));
  actionPBLinkSelection->setIconText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Link with selection", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
  actionPBLinkSelection->setToolTip(QCoreApplication::translate("pqPipelineBrowserContextMenu",
    "Link this source and current selected source as a selection link", Q_NULLPTR));
#endif // QT_NO_TOOLTIP

  QAction* actionPBRename = new QAction(menu.parent());
  actionPBRename->setObjectName(QStringLiteral("actionPBRename"));
  actionPBRename->setText(
    QCoreApplication::translate("pqPipelineBrowserContextMenu", "Rename", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
  actionPBRename->setToolTip(QCoreApplication::translate(
    "pqPipelineBrowserContextMenu", "Rename currently selected source", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_STATUSTIP
  actionPBRename->setStatusTip(QCoreApplication::translate(
    "pqPipelineBrowserContextMenu", "Rename currently selected source", Q_NULLPTR));
#endif // QT_NO_TOOLTIP

  menu.addAction(actionPBOpen);
  menu.addAction(actionPBShowAll);
  menu.addAction(actionPBHideAll);
  menu.addSeparator();
  menu.addAction(actionPBCopy);
  menu.addAction(actionPBPaste);
  menu.addAction(actionPBCopyPipeline);
  menu.addAction(actionPBPastePipeline);
  menu.addSeparator();
  menu.addAction(actionPBDelete);
  menu.addAction(actionPBDeleteTree);
  menu.addAction(actionPBRename);
  menu.addAction(actionPBReloadFiles);
  menu.addAction(actionPBChangeFile);
  menu.addAction(actionPBIgnoreTime);
  menu.addSeparator();
  menu.addAction(actionPBChangeInput);
  QMenu* addFilterMenu =
    menu.addMenu(QCoreApplication::translate("pqPipelineBrowserContextMenu", "Add Filter"));
  pqParaViewMenuBuilders::buildFiltersMenu(
    *addFilterMenu, nullptr, true /*hide disabled*/, false /*quickLaunchable*/);
  menu.addSeparator();
  menu.addAction(actionPBCreateCustomFilter);
  menu.addAction(actionPBLinkSelection);

  // And here the reactions come in handy! Just reuse the reaction used for
  // File | Open.
  new pqLoadDataReaction(actionPBOpen);
  new pqShowHideAllReaction(actionPBShowAll, pqShowHideAllReaction::ActionType::Show);
  new pqShowHideAllReaction(actionPBHideAll, pqShowHideAllReaction::ActionType::Hide);
  new pqCopyReaction(actionPBCopy);
  new pqCopyReaction(actionPBPaste, true);
  new pqCopyReaction(actionPBCopyPipeline, false, true);
  new pqCopyReaction(actionPBPastePipeline, true, true);
  new pqChangePipelineInputReaction(actionPBChangeInput);
  new pqReloadFilesReaction(actionPBReloadFiles);
  new pqIgnoreSourceTimeReaction(actionPBIgnoreTime);
  new pqDeleteReaction(actionPBDelete);
  new pqDeleteReaction(actionPBDeleteTree, pqDeleteReaction::DeleteModes::TREE);
  new pqCreateCustomFilterReaction(actionPBCreateCustomFilter);
  new pqLinkSelectionReaction(actionPBLinkSelection);
  new pqRenameProxyReaction(actionPBRename, mainWindow);
  new pqChangeFileNameReaction(actionPBChangeFile);
}

//-----------------------------------------------------------------------------
void pqParaViewMenuBuilders::buildMacrosMenu(QMenu& menu)
{
  Q_UNUSED(menu);
#if VTK_MODULE_ENABLE_ParaView_pqPython
  // Give the macros menu to the pqPythonMacroSupervisor
  pqPythonManager* manager = pqPVApplicationCore::instance()->pythonManager();
  if (manager)
  {
    new pqMacroReaction(
      menu.addAction(QCoreApplication::translate("pqMacrosMenu", "Import new macro..."))
      << pqSetName("actionMacroCreate"));

    new pqEditMacrosReaction(
      menu.addAction(QCoreApplication::translate("pqMacrosMenu", "Edit Macros"))
      << pqSetName("actionMacroEdit"));

    menu.addSeparator();
    manager->addWidgetForRunMacros(&menu);
  }

#endif
}

//-----------------------------------------------------------------------------
void pqParaViewMenuBuilders::buildHelpMenu(QMenu& menu)
{
  QString documentationPath(vtkPVFileInformation::GetParaViewDocDirectory().c_str());
  QString paraViewGettingStartedFile = documentationPath + "/GettingStarted.pdf";

  // Getting Started with ParaView
  new pqDesktopServicesReaction(QUrl::fromLocalFile(paraViewGettingStartedFile),
    (menu.addAction(QIcon(":/pqWidgets/Icons/pdf.png"),
       QCoreApplication::translate("pqHelpMenu", "Getting Started with ParaView"))
      << pqSetName("actionGettingStarted")));

  QString versionString = QString("%1.%2.%3")
                            .arg(vtkSMProxyManager::GetVersionMajor())
                            .arg(vtkSMProxyManager::GetVersionMinor())
                            .arg(vtkSMProxyManager::GetVersionPatch());

  // ParaView Guide. If there is a copy local to the install tree, use it instead of going
  // out to the web.
  QString paraViewGuideFile =
    QString("%1/ParaViewGuide-%2.pdf").arg(documentationPath).arg(versionString);
  QFile guideLocalFile(paraViewGuideFile);
  QAction* guide = menu.addAction(QCoreApplication::translate("pqHelpMenu", "ParaView Guide"));
  guide->setObjectName("actionGuide");
  guide->setShortcut(QKeySequence::HelpContents);
  if (guideLocalFile.exists())
  {
    guide->setIcon(QIcon(":/pqWidgets/Icons/pdf.png"));
    QUrl guideUrl = QUrl::fromLocalFile(paraViewGuideFile);
    new pqDesktopServicesReaction(guideUrl, guide);
  }
  else
  {
    // Remote ParaView Guide
    QString guideURL = QString("https://docs.paraview.org/en/v%1/").arg(versionString);
    new pqDesktopServicesReaction(QUrl(guideURL), guide);
  }

#ifdef PARAVIEW_USE_QTHELP
  // Help
  QAction* help = menu.addAction(QCoreApplication::translate(
                    "pqHelpMenu", "Reader, Filter, and Writer Reference"))
    << pqSetName("actionHelp");
  new pqHelpReaction(help);
#endif

  // -----------------
  menu.addSeparator();

  // ParaView Tutorial
  QString selfDirectedTutorialURL =
    QString("https://docs.paraview.org/en/v%1/Tutorials/SelfDirectedTutorial/index.html")
      .arg(versionString);
  new pqDesktopServicesReaction(QUrl(selfDirectedTutorialURL),
    (menu.addAction(QCoreApplication::translate("pqHelpMenu", "ParaView Self-directed Tutorial"))
      << pqSetName("actionSelfDirectedTutorialNotes")));

  // Classroom Tutorials by Sandia National Laboratories
  QString classroomTutorialsURL =
    QString("https://docs.paraview.org/en/v%1/Tutorials/ClassroomTutorials/index.html")
      .arg(versionString);
  new pqDesktopServicesReaction(QUrl(classroomTutorialsURL),
    (menu.addAction(QCoreApplication::translate("pqHelpMenu", "ParaView Classroom Tutorials"))
      << pqSetName("actionClassroomTutorials")));

  // Example Data Sets

  // Example Visualizations
  new pqExampleVisualizationsDialogReaction(
    menu.addAction(QCoreApplication::translate("pqHelpMenu", "Example Visualizations"))
    << pqSetName("ExampleVisualizations"));

  // -----------------
  menu.addSeparator();

  // ParaView Web Site
  new pqDesktopServicesReaction(QUrl("http://www.paraview.org"),
    (menu.addAction(QCoreApplication::translate("pqHelpMenu", "ParaView Web Site"))
      << pqSetName("actionWebSite")));

  // ParaView Wiki
  new pqDesktopServicesReaction(QUrl("http://www.paraview.org/Wiki/ParaView"),
    (menu.addAction(QCoreApplication::translate("pqHelpMenu", "ParaView Wiki"))
      << pqSetName("actionWiki")));

  // ParaView Community Support
  new pqDesktopServicesReaction(QUrl("https://discourse.paraview.org/c/paraview-support"),
    (menu.addAction(QCoreApplication::translate("pqHelpMenu", "ParaView Community Support"))
      << pqSetName("actionCommunitySupport")));

  // ParaView Release Notes
  versionString.replace('.', '-');
  new pqDesktopServicesReaction(
    QUrl("https://www.kitware.com/paraview-" + versionString + "-release-notes/"),
    (menu.addAction(QCoreApplication::translate("pqHelpMenu", "Release Notes"))
      << pqSetName("actionReleaseNotes")));

  // -----------------
  menu.addSeparator();

  // Professional Support
  new pqDesktopServicesReaction(QUrl("https://www.paraview.org/custom-support/"),
    (menu.addAction(QCoreApplication::translate("pqHelpMenu", "Professional Support"))
      << pqSetName("actionProfessionalSupport")));

  // Professional Training
  new pqDesktopServicesReaction(QUrl("https://www.kitware.com/support/#training"),
    (menu.addAction(QCoreApplication::translate("pqHelpMenu", "Professional Training"))
      << pqSetName("actionTraining")));

  // Online Tutorials
  new pqDesktopServicesReaction(QUrl("http://www.paraview.org/tutorials/"),
    (menu.addAction(QCoreApplication::translate("pqHelpMenu", "Online Tutorials"))
      << pqSetName("actionTutorials")));

  // Online Blogs
  new pqDesktopServicesReaction(QUrl("https://www.kitware.com/tag/ParaView/"),
    (menu.addAction(QCoreApplication::translate("pqHelpMenu", "Online Blogs"))
      << pqSetName("actionBlogs")));

#if !defined(__APPLE__)
  // -----------------
  menu.addSeparator();
#endif

  // Bug report
  QString bugReportURL = QString("https://gitlab.kitware.com/paraview/paraview/-/issues/"
                                 "new?issue[title]=Bug&issue[description]="
                                 "OS: %1\nParaView Version: %2\nQt Version: %3\n\nDescription")
                           .arg(QSysInfo::prettyProductName())
                           .arg(PARAVIEW_VERSION_FULL)
                           .arg(QT_VERSION_STR);
  new pqDesktopServicesReaction(QUrl(bugReportURL),
    (menu.addAction(QCoreApplication::translate("pqHelpMenu", "Bug Report"))
      << pqSetName("bugReport")));

  // About
  new pqAboutDialogReaction(menu.addAction(QCoreApplication::translate("pqHelpMenu", "About..."))
    << pqSetName("actionAbout"));
}

//-----------------------------------------------------------------------------
void pqParaViewMenuBuilders::buildToolbars(QMainWindow& mainWindow)
{
  QToolBar* mainToolBar = new pqMainControlsToolbar(&mainWindow)
    << pqSetName("MainControlsToolbar");
  mainToolBar->layout()->setSpacing(0);
  mainWindow.addToolBar(Qt::TopToolBarArea, mainToolBar);

  QToolBar* vcrToolbar = new pqVCRToolbar(&mainWindow) << pqSetName("VCRToolbar");
  vcrToolbar->layout()->setSpacing(0);
  mainWindow.addToolBar(Qt::TopToolBarArea, vcrToolbar);

  QToolBar* timeToolbar = new pqAnimationTimeToolbar(&mainWindow)
    << pqSetName("currentTimeToolbar");
  timeToolbar->layout()->setSpacing(0);
  mainWindow.addToolBar(Qt::TopToolBarArea, timeToolbar);

  auto* toolbarController = new pqCustomViewpointsDefaultController(&mainWindow);
  QToolBar* customViewpointsToolbar = new pqCustomViewpointsToolbar(toolbarController, &mainWindow)
    << pqSetName("customViewpointsToolbar");
  customViewpointsToolbar->layout()->setSpacing(0);
  mainWindow.addToolBar(Qt::TopToolBarArea, customViewpointsToolbar);

  mainWindow.addToolBarBreak();

  QToolBar* colorToolbar = new pqColorToolbar(&mainWindow) << pqSetName("variableToolbar");
  colorToolbar->layout()->setSpacing(0);
  mainWindow.addToolBar(Qt::TopToolBarArea, colorToolbar);

  QToolBar* reprToolbar = new pqRepresentationToolbar(&mainWindow)
    << pqSetName("representationToolbar");
  reprToolbar->layout()->setSpacing(0);
  mainWindow.addToolBar(Qt::TopToolBarArea, reprToolbar);

  QToolBar* cameraToolbar = new pqCameraToolbar(&mainWindow) << pqSetName("cameraToolbar");
  cameraToolbar->layout()->setSpacing(0);
  mainWindow.addToolBar(Qt::TopToolBarArea, cameraToolbar);

  QToolBar* axesToolbar = new pqAxesToolbar(&mainWindow) << pqSetName("axesToolbar");
  axesToolbar->layout()->setSpacing(0);
  mainWindow.addToolBar(Qt::TopToolBarArea, axesToolbar);

  QToolBar* lightingToolbar = new pqLightToolbar(&mainWindow) << pqSetName("lightingToolbar");
  lightingToolbar->layout()->setSpacing(0);
  mainWindow.addToolBar(Qt::TopToolBarArea, lightingToolbar);

#if VTK_MODULE_ENABLE_ParaView_pqPython
  // Give the macros menu to the pqPythonMacroSupervisor
  pqPythonManager* manager =
    qobject_cast<pqPythonManager*>(pqApplicationCore::instance()->manager("PYTHON_MANAGER"));
  if (manager)
  {
    QToolBar* macrosToolbar =
      new QToolBar(QCoreApplication::translate("pqmacrosToolbar", "Macros Toolbars"), &mainWindow)
      << pqSetName("MacrosToolbar");
    manager->addWidgetForRunMacros(macrosToolbar);
    mainWindow.addToolBar(Qt::TopToolBarArea, macrosToolbar);
  }
#endif
}

//-----------------------------------------------------------------------------
void pqParaViewMenuBuilders::buildCatalystMenu(QMenu& menu)
{
  new pqCatalystConnectReaction(
    menu.addAction(QCoreApplication::translate("pqCatalystMenu", "Connect..."))
    << pqSetName("actionCatalystConnect"));
  new pqCatalystPauseSimulationReaction(
    menu.addAction(QCoreApplication::translate("pqCatalystMenu", "Pause Simulation"))
    << pqSetName("actionCatalystPauseSimulation"));

  new pqCatalystContinueReaction(
    menu.addAction(QCoreApplication::translate("pqCatalystMenu", "Continue"))
    << pqSetName("actionCatalystContinue"));

  new pqCatalystSetBreakpointReaction(
    menu.addAction(QCoreApplication::translate("pqCatalystMenu", "Set Breakpoint"))
    << pqSetName("actionCatalystSetBreakpoint"));

  new pqCatalystRemoveBreakpointReaction(
    menu.addAction(QCoreApplication::translate("pqCatalystMenu", "Remove Breakpoint"))
    << pqSetName("actionCatalystRemoveBreakpoint"));
}
