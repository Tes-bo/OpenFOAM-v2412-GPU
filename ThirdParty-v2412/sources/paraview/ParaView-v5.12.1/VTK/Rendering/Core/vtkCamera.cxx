// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
#include "vtkCamera.h"

#include "vtkCallbackCommand.h"
#include "vtkInformation.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPerspectiveTransform.h"
#include "vtkRenderer.h"
#include "vtkTimeStamp.h"
#include "vtkTransform.h"

#include <cassert>
#include <cmath>

//------------------------------------------------------------------------------
// Needed when we don't use the vtkStandardNewMacro.
VTK_ABI_NAMESPACE_BEGIN
vtkObjectFactoryNewMacro(vtkCamera);

vtkCxxSetObjectMacro(vtkCamera, Information, vtkInformation);
vtkCxxSetObjectMacro(vtkCamera, EyeTransformMatrix, vtkMatrix4x4);
vtkCxxSetObjectMacro(vtkCamera, ModelTransformMatrix, vtkMatrix4x4);
vtkCxxSetObjectMacro(vtkCamera, ExplicitProjectionTransformMatrix, vtkMatrix4x4);

//------------------------------------------------------------------------------
class vtkCameraCallbackCommand : public vtkCommand
{
public:
  static vtkCameraCallbackCommand* New() { return new vtkCameraCallbackCommand; }
  vtkCamera* Self;
  void Execute(vtkObject*, unsigned long, void*) override
  {
    if (this->Self)
    {
      this->Self->Modified();
      this->Self->ComputeViewTransform();
      this->Self->ComputeDistance();
      this->Self->ComputeCameraLightTransform();
    }
  }

protected:
  vtkCameraCallbackCommand() { this->Self = nullptr; }
  ~vtkCameraCallbackCommand() override = default;
};

//------------------------------------------------------------------------------
// Construct camera instance with its focal point at the origin,
// and position=(0,0,1). The view up is along the y-axis,
// view angle is 30 degrees, and the clipping range is (.1,1000).
vtkCamera::vtkCamera()
{
  this->FocalPoint[0] = 0.0;
  this->FocalPoint[1] = 0.0;
  this->FocalPoint[2] = 0.0;

  this->Position[0] = 0.0;
  this->Position[1] = 0.0;
  this->Position[2] = 1.0;

  this->ViewUp[0] = 0.0;
  this->ViewUp[1] = 1.0;
  this->ViewUp[2] = 0.0;

  this->DirectionOfProjection[0] = 0.0;
  this->DirectionOfProjection[1] = 0.0;
  this->DirectionOfProjection[2] = 0.0;

  this->ViewAngle = 30.0;
  this->UseHorizontalViewAngle = 0;

  this->UseOffAxisProjection = 0;

  this->ScreenBottomLeft[0] = -0.5;
  this->ScreenBottomLeft[1] = -0.5;
  this->ScreenBottomLeft[2] = -0.5;

  this->ScreenBottomRight[0] = 0.5;
  this->ScreenBottomRight[1] = -0.5;
  this->ScreenBottomRight[2] = -0.5;

  this->ScreenTopRight[0] = 0.5;
  this->ScreenTopRight[1] = 0.5;
  this->ScreenTopRight[2] = -0.5;

  this->ScreenCenter[0] = 0.0;
  this->ScreenCenter[1] = 0.0;
  this->ScreenCenter[2] = -0.5;

  this->OffAxisClippingAdjustment = 0.0;
  this->EyeSeparation = 0.06;

  this->EyeTransformMatrix = vtkMatrix4x4::New();
  this->EyeTransformMatrix->Identity();

  this->ProjectionPlaneOrientationMatrix = nullptr;

  this->ModelTransformMatrix = vtkMatrix4x4::New();
  this->ModelTransformMatrix->Identity();

  this->ClippingRange[0] = 0.01;
  this->ClippingRange[1] = 1000.01;
  this->Thickness = 1000.0;

  this->ParallelProjection = 0;
  this->ParallelScale = 1.0;

  this->EyeAngle = 2.0;
  this->Stereo = 0;
  this->LeftEye = 1;

  this->WindowCenter[0] = 0.0;
  this->WindowCenter[1] = 0.0;

  this->ViewShear[0] = 0.0;
  this->ViewShear[1] = 0.0;
  this->ViewShear[2] = 1.0;

  this->FocalDisk = 1.0;
  this->FocalDistance = 0.0;

  this->Transform = vtkPerspectiveTransform::New();
  this->ViewTransform = vtkTransform::New();
  this->ProjectionTransform = vtkPerspectiveTransform::New();
  this->CameraLightTransform = vtkTransform::New();
  this->ModelViewTransform = vtkTransform::New();
  this->ExplicitProjectionTransformMatrix = nullptr;
  this->UseExplicitProjectionTransformMatrix = false;
  this->UserTransform = nullptr;
  this->UserViewTransform = nullptr;
  this->UserViewTransformCallbackCommand = nullptr;

  // initialize the ViewTransform
  this->ComputeViewTransform();
  this->ComputeDistance();
  this->ComputeCameraLightTransform();

  this->FreezeFocalPoint = false;
  this->UseScissor = false;

  this->Information = vtkInformation::New();
  this->Information->Register(this);
  this->Information->Delete();

  this->ExplicitAspectRatio = 1.0;
  this->UseExplicitAspectRatio = false;

  this->FocalPointScale = 1.0;
  this->FocalPointShift[0] = 0.0;
  this->FocalPointShift[1] = 0.0;
  this->FocalPointShift[2] = 0.0;
  this->NearPlaneScale = 1.0;
  this->NearPlaneShift[0] = 0.0;
  this->NearPlaneShift[1] = 0.0;
  this->NearPlaneShift[2] = 0.0;
  this->ShiftScaleThreshold = 2.0;
}

//------------------------------------------------------------------------------
vtkCamera::~vtkCamera()
{
  this->EyeTransformMatrix->Delete();
  this->EyeTransformMatrix = nullptr;

  if (this->ProjectionPlaneOrientationMatrix)
  {
    this->ProjectionPlaneOrientationMatrix->Delete();
    this->ProjectionPlaneOrientationMatrix = nullptr;
  }

  this->ModelTransformMatrix->Delete();
  this->ModelTransformMatrix = nullptr;

  this->Transform->Delete();
  this->ViewTransform->Delete();
  this->ProjectionTransform->Delete();
  this->CameraLightTransform->Delete();
  this->ModelViewTransform->Delete();
  if (this->ExplicitProjectionTransformMatrix)
  {
    this->ExplicitProjectionTransformMatrix->UnRegister(this);
    this->ExplicitProjectionTransformMatrix = nullptr;
  }
  if (this->UserTransform)
  {
    this->UserTransform->UnRegister(this);
    this->UserTransform = nullptr;
  }
  if (this->UserViewTransform)
  {
    this->UserViewTransform->RemoveObserver(this->UserViewTransformCallbackCommand);
    this->UserViewTransform->UnRegister(this);
    this->UserViewTransform = nullptr;
  }
  if (this->UserViewTransformCallbackCommand)
  {
    this->UserViewTransformCallbackCommand->Delete();
  }

  this->SetInformation(nullptr);
}

//------------------------------------------------------------------------------
void vtkCamera::SetScissorRect(vtkRecti scissorRect)
{
  this->ScissorRect = scissorRect;
}

//------------------------------------------------------------------------------
void vtkCamera::GetScissorRect(vtkRecti& scissorRect)
{
  scissorRect = this->ScissorRect;
}

//------------------------------------------------------------------------------
// The first set of methods deal exclusively with the ViewTransform, which
// is the only transform which is set up entirely in the camera.  The
// perspective transform must be set up by the Renderer because the
// Camera doesn't know the Renderer's aspect ratio.
//------------------------------------------------------------------------------
void vtkCamera::SetPosition(double x, double y, double z)
{
  if (x == this->Position[0] && y == this->Position[1] && z == this->Position[2])
  {
    return;
  }

  this->Position[0] = x;
  this->Position[1] = y;
  this->Position[2] = z;

  vtkDebugMacro(<< " Position set to ( " << this->Position[0] << ", " << this->Position[1] << ", "
                << this->Position[2] << ")");

  this->ComputeViewTransform();
  // recompute the focal distance
  this->ComputeDistance();
  this->ComputeCameraLightTransform();

  this->Modified();
}

//------------------------------------------------------------------------------
void vtkCamera::SetUserTransform(vtkHomogeneousTransform* transform)
{
  if (transform == this->UserTransform)
  {
    return;
  }
  if (this->UserTransform)
  {
    this->UserTransform->Delete();
    this->UserTransform = nullptr;
  }
  if (transform)
  {
    this->UserTransform = transform;
    this->UserTransform->Register(this);
  }
  this->Modified();
}

//------------------------------------------------------------------------------
void vtkCamera::SetUserViewTransform(vtkHomogeneousTransform* transform)
{
  if (transform == this->UserViewTransform)
  {
    return;
  }
  if (this->UserViewTransform)
  {
    this->UserViewTransform->RemoveObserver(this->UserViewTransformCallbackCommand);
    this->UserViewTransform->Delete();
    this->UserViewTransform = nullptr;
  }
  if (transform)
  {
    this->UserViewTransform = transform;
    this->UserViewTransform->Register(this);
    if (!this->UserViewTransformCallbackCommand)
    {
      this->UserViewTransformCallbackCommand = vtkCameraCallbackCommand::New();
      this->UserViewTransformCallbackCommand->Self = this;
    }
    this->UserViewTransform->AddObserver(
      vtkCommand::ModifiedEvent, this->UserViewTransformCallbackCommand);
  }
  this->Modified();
  this->ComputeViewTransform();
  this->ComputeDistance();
  this->ComputeCameraLightTransform();
}

//------------------------------------------------------------------------------
void vtkCamera::SetFocalPoint(double x, double y, double z)
{
  if (x == this->FocalPoint[0] && y == this->FocalPoint[1] && z == this->FocalPoint[2])
  {
    return;
  }

  this->FocalPoint[0] = x;
  this->FocalPoint[1] = y;
  this->FocalPoint[2] = z;

  vtkDebugMacro(<< " FocalPoint set to ( " << this->FocalPoint[0] << ", " << this->FocalPoint[1]
                << ", " << this->FocalPoint[2] << ")");

  this->ComputeViewTransform();
  // recompute the focal distance
  this->ComputeDistance();
  this->ComputeCameraLightTransform();

  this->Modified();
}

//------------------------------------------------------------------------------
void vtkCamera::SetViewUp(double x, double y, double z)
{
  // normalize ViewUp, but do _not_ orthogonalize it by default
  double norm = sqrt(x * x + y * y + z * z);

  if (norm != 0)
  {
    x /= norm;
    y /= norm;
    z /= norm;
  }
  else
  {
    x = 0;
    y = 1;
    z = 0;
  }

  if (x == this->ViewUp[0] && y == this->ViewUp[1] && z == this->ViewUp[2])
  {
    return;
  }

  this->ViewUp[0] = x;
  this->ViewUp[1] = y;
  this->ViewUp[2] = z;

  vtkDebugMacro(<< " ViewUp set to ( " << this->ViewUp[0] << ", " << this->ViewUp[1] << ", "
                << this->ViewUp[2] << ")");

  this->ComputeViewTransform();
  this->ComputeCameraLightTransform();
  this->Modified();
}

//------------------------------------------------------------------------------
// The ViewTransform depends on only three ivars:  the Position, the
// FocalPoint, and the ViewUp vector.  All the other methods are there
// simply for the sake of the users' convenience.
void vtkCamera::ComputeViewTransform()
{
  // main view through the camera
  this->Transform->Identity();
  if (this->UserViewTransform)
  {
    this->Transform->Concatenate(this->UserViewTransform);
  }
  if (!this->UseOffAxisProjection)
  {
    this->Transform->SetupCamera(this->Position, this->FocalPoint, this->ViewUp);
  }
  else
  {
    double pe[3] = { 0.0 };
    this->GetStereoEyePosition(pe);

    // Create the view point offset matrix
    vtkSmartPointer<vtkMatrix4x4> T = vtkSmartPointer<vtkMatrix4x4>::New();
    T->SetElement(0, 3, -pe[0]);
    T->SetElement(1, 3, -pe[1]);
    T->SetElement(2, 3, -pe[2]);

    this->Transform->Concatenate(T);
  }
  this->ViewTransform->Identity();
  this->ViewTransform->Concatenate(this->Transform->GetMatrix());
}

//------------------------------------------------------------------------------
void vtkCamera::ComputeCameraLightTransform()
{
  vtkTransform* t;
  double d;

  // assumes a valid view transform and valid camera distance

  t = this->CameraLightTransform;
  t->Identity();
  t->SetMatrix(this->ViewTransform->GetMatrix());
  t->Inverse();

  d = this->Distance;
  t->Scale(d, d, d);
  t->Translate(0.0, 0.0, -1.0);
}

//------------------------------------------------------------------------------
void vtkCamera::ComputeScreenOrientationMatrix()
{
  if (this->ProjectionPlaneOrientationMatrix == nullptr)
  {
    // Compute the screen orientation matrix lazily, and only once.
    double vr[3] = { 0.0 };
    double vu[3] = { 0.0 };
    double vn[3] = { 0.0 };
    double screenDiag[3] = { 0 };

    for (int i = 0; i < 3; ++i)
    {
      vr[i] = this->ScreenBottomRight[i] - this->ScreenBottomLeft[i];
      vu[i] = this->ScreenTopRight[i] - this->ScreenBottomRight[i];

      this->ScreenCenter[i] = (this->ScreenBottomLeft[i] + this->ScreenTopRight[i]) / 2.0;
      screenDiag[i] = this->ScreenBottomLeft[i] - this->ScreenTopRight[i];
    }

    this->OffAxisClippingAdjustment = vtkMath::Norm(screenDiag);

    vtkMath::Normalize(vr);
    vtkMath::Normalize(vu);
    vtkMath::Cross(vr, vu, vn);
    vtkMath::Normalize(vn);

    this->ProjectionPlaneOrientationMatrix = vtkMatrix4x4::New();

    this->ProjectionPlaneOrientationMatrix->SetElement(0, 0, vr[0]);
    this->ProjectionPlaneOrientationMatrix->SetElement(0, 1, vr[1]);
    this->ProjectionPlaneOrientationMatrix->SetElement(0, 2, vr[2]);
    this->ProjectionPlaneOrientationMatrix->SetElement(1, 0, vu[0]);
    this->ProjectionPlaneOrientationMatrix->SetElement(1, 1, vu[1]);
    this->ProjectionPlaneOrientationMatrix->SetElement(1, 2, vu[2]);
    this->ProjectionPlaneOrientationMatrix->SetElement(2, 0, vn[0]);
    this->ProjectionPlaneOrientationMatrix->SetElement(2, 1, vn[1]);
    this->ProjectionPlaneOrientationMatrix->SetElement(2, 2, vn[2]);
  }
}

//------------------------------------------------------------------------------
void vtkCamera::ComputeOffAxisProjectionFrustum()
{
  // This version of off-axis projection was implemented from the article
  // referenced below, and variable names in this method were chosen to match
  // those used in the article.
  //
  // TItle: Generalized perspective projection
  // Author: Robert Kooima
  // Date: 2009/6
  // Journal: J. Sch. Electron. Eng. Comput. Sci
  // Volume: 6

  this->ComputeScreenOrientationMatrix();

  double n = this->ClippingRange[0];
  double f = this->ClippingRange[1];
  double pe[3] = { 0.0 };
  this->GetStereoEyePosition(pe);

  double pa[4] = { this->ScreenBottomLeft[0], this->ScreenBottomLeft[1], this->ScreenBottomLeft[2],
    1.0 };
  double pb[4] = { this->ScreenBottomRight[0], this->ScreenBottomRight[1],
    this->ScreenBottomRight[2], 1.0 };
  double pc[4] = { this->ScreenTopRight[0], this->ScreenTopRight[1], this->ScreenTopRight[2], 1.0 };

  double vr[3];
  vr[0] = this->ProjectionPlaneOrientationMatrix->GetElement(0, 0);
  vr[1] = this->ProjectionPlaneOrientationMatrix->GetElement(0, 1);
  vr[2] = this->ProjectionPlaneOrientationMatrix->GetElement(0, 2);

  double vu[3];
  vu[0] = this->ProjectionPlaneOrientationMatrix->GetElement(1, 0);
  vu[1] = this->ProjectionPlaneOrientationMatrix->GetElement(1, 1);
  vu[2] = this->ProjectionPlaneOrientationMatrix->GetElement(1, 2);

  double vn[3] = { 0.0 };
  vn[0] = this->ProjectionPlaneOrientationMatrix->GetElement(2, 0);
  vn[1] = this->ProjectionPlaneOrientationMatrix->GetElement(2, 1);
  vn[2] = this->ProjectionPlaneOrientationMatrix->GetElement(2, 2);

  double va[3] = { 0.0 };
  double vb[3] = { 0.0 };
  double vc[3] = { 0.0 };

  for (int i = 0; i < 3; ++i)
  {
    va[i] = pa[i] - pe[i];
    vb[i] = pb[i] - pe[i];
    vc[i] = pc[i] - pe[i];
  }

  double d = -vtkMath::Dot(vn, va);
  double nOverD = n / d;

  double l = vtkMath::Dot(vr, va) * (nOverD);
  double r = vtkMath::Dot(vr, vb) * (nOverD);
  double b = vtkMath::Dot(vu, va) * (nOverD);
  double t = vtkMath::Dot(vu, vc) * (nOverD);

  // Populate it as glFrustum would do
  this->ProjectionTransform->GetMatrix()->SetElement(0, 0, (2.0 * n) / (r - l));
  this->ProjectionTransform->GetMatrix()->SetElement(0, 1, 0.0);
  this->ProjectionTransform->GetMatrix()->SetElement(0, 2, (r + l) / (r - l));
  this->ProjectionTransform->GetMatrix()->SetElement(0, 3, 0.0);

  this->ProjectionTransform->GetMatrix()->SetElement(1, 0, 0.0);
  this->ProjectionTransform->GetMatrix()->SetElement(1, 1, (2.0 * n) / (t - b));
  this->ProjectionTransform->GetMatrix()->SetElement(1, 2, (t + b) / (t - b));
  this->ProjectionTransform->GetMatrix()->SetElement(1, 3, 0.0);

  this->ProjectionTransform->GetMatrix()->SetElement(2, 0, 0.0);
  this->ProjectionTransform->GetMatrix()->SetElement(2, 1, 0.0);
  this->ProjectionTransform->GetMatrix()->SetElement(2, 2, -(f + n) / (f - n));
  this->ProjectionTransform->GetMatrix()->SetElement(2, 3, -(2.0 * f * n) / (f - n));

  this->ProjectionTransform->GetMatrix()->SetElement(3, 0, 0.0);
  this->ProjectionTransform->GetMatrix()->SetElement(3, 1, 0.0);
  this->ProjectionTransform->GetMatrix()->SetElement(3, 2, -1.0);
  this->ProjectionTransform->GetMatrix()->SetElement(3, 3, 0.0);

  vtkMatrix4x4::Multiply4x4(this->ProjectionTransform->GetMatrix(),
    this->ProjectionPlaneOrientationMatrix, this->ProjectionTransform->GetMatrix());

  // The viewer offset translation matrix, T, described in the paper, is kept
  // in the view transform (see ComputeViewTransform()).  It's important to keep
  // it there for lighting purposes.
}

//------------------------------------------------------------------------------
double vtkCamera::GetOffAxisClippingAdjustment()
{
  return this->OffAxisClippingAdjustment;
}

//------------------------------------------------------------------------------
void vtkCamera::ComputeModelViewMatrix()
{
  if (this->ModelViewTransform->GetMTime() < this->ModelTransformMatrix->GetMTime() ||
    this->ModelViewTransform->GetMTime() < this->ViewTransform->GetMTime() ||
    (this->UseOffAxisProjection &&
      this->ModelViewTransform->GetMTime() < this->EyeTransformMatrix->GetMTime()))
  {
    if (this->UseOffAxisProjection)
    {
      this->ComputeViewTransform();
    }
    vtkMatrix4x4::Multiply4x4(this->ViewTransform->GetMatrix(), this->ModelTransformMatrix,
      this->ModelViewTransform->GetMatrix());
  }
}

//------------------------------------------------------------------------------
void vtkCamera::OrthogonalizeViewUp()
{
  // the orthogonalized ViewUp is just the second row of the view matrix
  vtkMatrix4x4* matrix = this->ViewTransform->GetMatrix();
  this->ViewUp[0] = matrix->GetElement(1, 0);
  this->ViewUp[1] = matrix->GetElement(1, 1);
  this->ViewUp[2] = matrix->GetElement(1, 2);

  this->Modified();
}

//------------------------------------------------------------------------------
// Set the distance of the focal point from the camera. The focal point is
// modified accordingly. This should be positive.
void vtkCamera::SetDistance(double d)
{
  if (this->Distance == d)
  {
    return;
  }

  this->Distance = d;

  // Distance should be greater than .0002
  if (this->Distance < 0.0002)
  {
    this->Distance = 0.0002;
    vtkDebugMacro(<< " Distance is set to minimum.");
  }

  // we want to keep the camera pointing in the same direction
  double* vec = this->DirectionOfProjection;

  // recalculate FocalPoint
  this->FocalPoint[0] = this->Position[0] + vec[0] * this->Distance;
  this->FocalPoint[1] = this->Position[1] + vec[1] * this->Distance;
  this->FocalPoint[2] = this->Position[2] + vec[2] * this->Distance;

  vtkDebugMacro(<< " Distance set to ( " << this->Distance << ")");

  this->ComputeViewTransform();
  this->ComputeCameraLightTransform();
  this->Modified();
}

//------------------------------------------------------------------------------
// This method must be called when the focal point or camera position changes
void vtkCamera::ComputeDistance()
{
  double dx = this->FocalPoint[0] - this->Position[0];
  double dy = this->FocalPoint[1] - this->Position[1];
  double dz = this->FocalPoint[2] - this->Position[2];

  this->Distance = sqrt(dx * dx + dy * dy + dz * dz);

  if (this->Distance < 1e-20)
  {
    this->Distance = 1e-20;
    vtkDebugMacro(<< " Distance is set to minimum.");

    double* vec = this->DirectionOfProjection;

    // recalculate FocalPoint
    this->FocalPoint[0] = this->Position[0] + vec[0] * this->Distance;
    this->FocalPoint[1] = this->Position[1] + vec[1] * this->Distance;
    this->FocalPoint[2] = this->Position[2] + vec[2] * this->Distance;
  }

  this->DirectionOfProjection[0] = dx / this->Distance;
  this->DirectionOfProjection[1] = dy / this->Distance;
  this->DirectionOfProjection[2] = dz / this->Distance;

  this->ComputeViewPlaneNormal();
}

//------------------------------------------------------------------------------
// Move the position of the camera along the view plane normal. Moving
// towards the focal point (e.g., > 1) is a dolly-in, moving away
// from the focal point (e.g., < 1) is a dolly-out.
void vtkCamera::Dolly(double amount)
{
  if (amount <= 0.0)
  {
    return;
  }

  // dolly moves the camera towards the focus
  double d = this->Distance / amount;

  this->SetPosition(this->FocalPoint[0] - d * this->DirectionOfProjection[0],
    this->FocalPoint[1] - d * this->DirectionOfProjection[1],
    this->FocalPoint[2] - d * this->DirectionOfProjection[2]);
}

//------------------------------------------------------------------------------
// Set the roll angle of the camera about the direction of projection
void vtkCamera::SetRoll(double roll)
{
  // roll is a rotation of camera view up about the direction of projection
  vtkDebugMacro(<< " Setting Roll to " << roll << "");

  // subtract the current roll
  roll -= this->GetRoll();

  if (fabs(roll) < 0.00001)
  {
    return;
  }

  this->Roll(roll);
}

//------------------------------------------------------------------------------
// Returns the roll of the camera.
double vtkCamera::GetRoll()
{
  double orientation[3];
  this->ViewTransform->GetOrientation(orientation);
  return orientation[2];
}

//------------------------------------------------------------------------------
// Rotate the camera around the view plane normal.
void vtkCamera::Roll(double angle)
{
  double newViewUp[3];
  this->Transform->Identity();

  // rotate ViewUp about the Direction of Projection
  this->Transform->RotateWXYZ(angle, this->DirectionOfProjection);

  // okay, okay, TransformPoint shouldn't be used on vectors -- but
  // the transform is rotation with no translation so this works fine.
  this->Transform->TransformPoint(this->ViewUp, newViewUp);
  this->SetViewUp(newViewUp);
}

//------------------------------------------------------------------------------
// Rotate the focal point about the view up vector centered at the camera's
// position.
void vtkCamera::Yaw(double angle)
{
  double newFocalPoint[3];
  double* pos = this->Position;
  this->Transform->Identity();

  // translate the camera to the origin,
  // rotate about axis,
  // translate back again
  this->Transform->Translate(+pos[0], +pos[1], +pos[2]);
  this->Transform->RotateWXYZ(angle, this->ViewUp);
  this->Transform->Translate(-pos[0], -pos[1], -pos[2]);

  // now transform focal point
  this->Transform->TransformPoint(this->FocalPoint, newFocalPoint);
  this->SetFocalPoint(newFocalPoint);
}

//------------------------------------------------------------------------------
// Rotate the focal point about the cross product of the view up vector
// and the negative of the , centered at the camera's position.
void vtkCamera::Pitch(double angle)
{
  double axis[3], newFocalPoint[3], savedViewUp[3];
  double* pos = this->Position;
  this->Transform->Identity();

  // the axis is the first row of the view transform matrix
  axis[0] = this->ViewTransform->GetMatrix()->GetElement(0, 0);
  axis[1] = this->ViewTransform->GetMatrix()->GetElement(0, 1);
  axis[2] = this->ViewTransform->GetMatrix()->GetElement(0, 2);

  // temporarily set the view up with the transformation applied
  // to avoid bad cross product computations during SetFocalPoint call
  this->GetViewUp(savedViewUp);
  this->Transform->RotateWXYZ(angle, axis);
  this->Transform->TransformPoint(this->ViewUp, this->ViewUp);
  this->Transform->Identity();

  // translate the camera to the origin,
  // rotate about axis,
  // translate back again
  this->Transform->Translate(+pos[0], +pos[1], +pos[2]);
  this->Transform->RotateWXYZ(angle, axis);
  this->Transform->Translate(-pos[0], -pos[1], -pos[2]);

  // now transform focal point
  this->Transform->TransformPoint(this->FocalPoint, newFocalPoint);
  this->SetFocalPoint(newFocalPoint);

  // restore the previous ViewUp vector
  this->ViewUp[0] = savedViewUp[0];
  this->ViewUp[1] = savedViewUp[1];
  this->ViewUp[2] = savedViewUp[2];
  // this is needed since the last time Modified was called (in SetFocalPoint),
  // the ViewUp was not same as savedViewUp. Since we're changing its value
  // here, we need to fire Modified event. We don't call `SetViewUp` since we
  // don't want the computation of the view transform to happen again.
  this->Modified();
}

//------------------------------------------------------------------------------
// Rotate the camera about the view up vector centered at the focal point.
void vtkCamera::Azimuth(double angle)
{
  double newPosition[3];
  double* fp = this->FocalPoint;
  this->Transform->Identity();

  // translate the focal point to the origin,
  // rotate about view up,
  // translate back again
  this->Transform->Translate(+fp[0], +fp[1], +fp[2]);
  this->Transform->RotateWXYZ(angle, this->ViewUp);
  this->Transform->Translate(-fp[0], -fp[1], -fp[2]);

  // apply the transform to the position
  this->Transform->TransformPoint(this->Position, newPosition);
  this->SetPosition(newPosition);
}

//------------------------------------------------------------------------------
// Rotate the camera about the cross product of the negative of the
// direction of projection and the view up vector centered on the focal point.
void vtkCamera::Elevation(double angle)
{
  double axis[3], newPosition[3], savedViewUp[3];
  double* fp = this->FocalPoint;
  this->Transform->Identity();

  // snatch the axis from the view transform matrix
  axis[0] = -this->ViewTransform->GetMatrix()->GetElement(0, 0);
  axis[1] = -this->ViewTransform->GetMatrix()->GetElement(0, 1);
  axis[2] = -this->ViewTransform->GetMatrix()->GetElement(0, 2);

  // temporarily set the view up with the transformation applied
  // to avoid bad cross product computations during SetPosition call
  this->GetViewUp(savedViewUp);
  this->Transform->RotateWXYZ(angle, axis);
  this->Transform->TransformPoint(this->ViewUp, this->ViewUp);
  this->Transform->Identity();

  // translate the focal point to the origin,
  // rotate about axis,
  // translate back again
  this->Transform->Translate(+fp[0], +fp[1], +fp[2]);
  this->Transform->RotateWXYZ(angle, axis);
  this->Transform->Translate(-fp[0], -fp[1], -fp[2]);

  // now transform position
  this->Transform->TransformPoint(this->Position, newPosition);
  this->SetPosition(newPosition);

  // restore the previous ViewUp vector
  this->ViewUp[0] = savedViewUp[0];
  this->ViewUp[1] = savedViewUp[1];
  this->ViewUp[2] = savedViewUp[2];
  // this is needed since the last time Modified was called (in SetPosition),
  // the ViewUp was not same as savedViewUp. Since we're changing its value
  // here, we need to fire Modified event. We don't call `SetViewUp` since we
  // don't want the computation of the view transform to happen again.
  this->Modified();
}

//------------------------------------------------------------------------------
// Apply Transform to camera
void vtkCamera::ApplyTransform(vtkTransform* t)
{
  double posOld[4], posNew[4], fpOld[4], fpNew[4], vuOld[4], vuNew[4];

  this->GetPosition(posOld);
  this->GetFocalPoint(fpOld);
  this->GetViewUp(vuOld);

  posOld[3] = 1.0;
  fpOld[3] = 1.0;
  vuOld[3] = 1.0;

  vuOld[0] += posOld[0];
  vuOld[1] += posOld[1];
  vuOld[2] += posOld[2];

  t->MultiplyPoint(posOld, posNew);
  t->MultiplyPoint(fpOld, fpNew);
  t->MultiplyPoint(vuOld, vuNew);

  vuNew[0] -= posNew[0];
  vuNew[1] -= posNew[1];
  vuNew[2] -= posNew[2];

  this->SetPosition(posNew);
  this->SetFocalPoint(fpNew);
  this->SetViewUp(vuNew);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// The following methods set up the information that the Renderer needs
// to set up the perspective transform.  The transformation matrix is
// created using the GetPerspectiveTransformMatrix method.
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void vtkCamera::SetParallelProjection(vtkTypeBool flag)
{
  if (this->ParallelProjection != flag)
  {
    this->ParallelProjection = flag;
    this->Modified();
    this->ViewingRaysModified();
  }
}

//------------------------------------------------------------------------------
void vtkCamera::SetViewAngle(double angle)
{
  double min = 0.00000001;
  double max = 179.0;

  if (this->ViewAngle != angle)
  {
    this->ViewAngle = (angle < min ? min : (angle > max ? max : angle));
    this->Modified();
    this->ViewingRaysModified();
  }
}

//------------------------------------------------------------------------------
void vtkCamera::SetUseHorizontalViewAngle(vtkTypeBool flag)
{
  if (flag == this->UseHorizontalViewAngle)
  {
    return;
  }
  this->UseHorizontalViewAngle = flag;
  this->Modified();
  this->ViewingRaysModified();
}

//------------------------------------------------------------------------------
void vtkCamera::SetParallelScale(double scale)
{
  if (this->ParallelScale != scale)
  {
    this->ParallelScale = scale;
    this->Modified();
    this->ViewingRaysModified();
  }
}

//------------------------------------------------------------------------------
// Change the ViewAngle (for perspective) or the ParallelScale (for parallel)
// so that more or less of a scene occupies the viewport.  A value > 1 is a
// zoom-in. A value < 1 is a zoom-out.
void vtkCamera::Zoom(double amount)
{
  if (amount <= 0.0)
  {
    return;
  }

  if (this->ParallelProjection)
  {
    this->SetParallelScale(this->ParallelScale / amount);
  }
  else
  {
    this->SetViewAngle(this->ViewAngle / amount);
  }
}

//------------------------------------------------------------------------------
void vtkCamera::SetClippingRange(double nearz, double farz)
{
  double thickness;

  // check the order
  if (nearz > farz)
  {
    vtkDebugMacro(<< " Front and back clipping range reversed");
    double temp = nearz;
    nearz = farz;
    farz = temp;
  }

  thickness = farz - nearz;

  // thickness should be greater than 1e-20
  if (thickness < 1e-20)
  {
    thickness = 1e-20;
    vtkDebugMacro(<< " ClippingRange thickness is set to minimum.");

    // set back plane
    farz = nearz + thickness;
  }

  if (nearz == this->ClippingRange[0] && farz == this->ClippingRange[1] &&
    this->Thickness == thickness)
  {
    return;
  }

  this->ClippingRange[0] = nearz;
  this->ClippingRange[1] = farz;
  this->Thickness = thickness;

  vtkDebugMacro(<< " ClippingRange set to ( " << this->ClippingRange[0] << ", "
                << this->ClippingRange[1] << ")");

  this->Modified();
}

//------------------------------------------------------------------------------
// Set the distance between clipping planes.
// This method adjusts the back clipping plane to the specified thickness
// behind the front clipping plane
void vtkCamera::SetThickness(double s)
{
  if (this->Thickness == s)
  {
    return;
  }

  this->Thickness = s;

  // thickness should be greater than 1e-20
  if (this->Thickness < 1e-20)
  {
    this->Thickness = 1e-20;
    vtkDebugMacro(<< " ClippingRange thickness is set to minimum.");
  }

  // set back plane
  this->ClippingRange[1] = this->ClippingRange[0] + this->Thickness;

  vtkDebugMacro(<< " ClippingRange set to ( " << this->ClippingRange[0] << ", "
                << this->ClippingRange[1] << ")");

  this->Modified();
}

//------------------------------------------------------------------------------
void vtkCamera::SetWindowCenter(double x, double y)
{
  if (this->WindowCenter[0] != x || this->WindowCenter[1] != y)
  {
    this->Modified();
    this->ViewingRaysModified();
    this->WindowCenter[0] = x;
    this->WindowCenter[1] = y;
  }
}

//------------------------------------------------------------------------------
void vtkCamera::SetObliqueAngles(double alpha, double beta)
{
  alpha = vtkMath::RadiansFromDegrees(alpha);
  beta = vtkMath::RadiansFromDegrees(beta);

  double cotbeta = cos(beta) / sin(beta);
  double dxdz = cos(alpha) * cotbeta;
  double dydz = sin(alpha) * cotbeta;

  this->SetViewShear(dxdz, dydz, 1.0);
}

//------------------------------------------------------------------------------
// Set the shear transform of the viewing frustum.  Parameters are
// dx/dz, dy/dz, and center.  center is a factor that describes where
// to shear around. The distance dshear from the camera where
// no shear occurs is given by (dshear = center * FocalDistance).
//
void vtkCamera::SetViewShear(double dxdz, double dydz, double center)
{
  if (dxdz != this->ViewShear[0] || dydz != this->ViewShear[1] || center != this->ViewShear[2])
  {
    this->Modified();
    this->ViewingRaysModified();

    this->ViewShear[0] = dxdz;
    this->ViewShear[1] = dydz;
    this->ViewShear[2] = center;

    this->ComputeViewPlaneNormal();
  }
}
//------------------------------------------------------------------------------

void vtkCamera::SetViewShear(double d[3])
{
  this->SetViewShear(d[0], d[1], d[2]);
}

//------------------------------------------------------------------------------
// Compute the projection transform matrix. This is used in converting
// between view and world coordinates.
void vtkCamera::ComputeProjectionTransform(double aspect, double nearz, double farz)
{
  this->ProjectionTransform->Identity();

  // apply user defined transform last if there is one
  if (this->UserTransform)
  {
    this->ProjectionTransform->Concatenate(this->UserTransform->GetMatrix());
  }

  if (this->UseExplicitProjectionTransformMatrix)
  {
    assert(this->ExplicitProjectionTransformMatrix != nullptr);
    this->ProjectionTransform->Concatenate(this->ExplicitProjectionTransformMatrix);
    return;
  }

  if (this->UseExplicitAspectRatio)
  {
    aspect = this->ExplicitAspectRatio;
  }

  // adjust Z-buffer range
  this->ProjectionTransform->AdjustZBuffer(-1, +1, nearz, farz);

  if (this->ParallelProjection)
  {
    // set up a rectangular parallelipiped

    double width = this->ParallelScale * aspect;
    double height = this->ParallelScale;

    double xmin = (this->WindowCenter[0] - 1.0) * width;
    double xmax = (this->WindowCenter[0] + 1.0) * width;
    double ymin = (this->WindowCenter[1] - 1.0) * height;
    double ymax = (this->WindowCenter[1] + 1.0) * height;

    this->ProjectionTransform->Ortho(
      xmin, xmax, ymin, ymax, this->ClippingRange[0], this->ClippingRange[1]);
  }
  else if (this->UseOffAxisProjection)
  {
    this->ComputeOffAxisProjectionFrustum();
  }
  else
  {
    // set up a perspective frustum

    double tmp = tan(vtkMath::RadiansFromDegrees(this->ViewAngle) / 2.);
    double width;
    double height;
    if (this->UseHorizontalViewAngle)
    {
      width = this->ClippingRange[0] * tmp;
      height = this->ClippingRange[0] * tmp / aspect;
    }
    else
    {
      width = this->ClippingRange[0] * tmp * aspect;
      height = this->ClippingRange[0] * tmp;
    }

    double xmin = (this->WindowCenter[0] - 1.0) * width;
    double xmax = (this->WindowCenter[0] + 1.0) * width;
    double ymin = (this->WindowCenter[1] - 1.0) * height;
    double ymax = (this->WindowCenter[1] + 1.0) * height;

    this->ProjectionTransform->Frustum(
      xmin, xmax, ymin, ymax, this->ClippingRange[0], this->ClippingRange[1]);
  }

  if (this->Stereo && !this->UseOffAxisProjection)
  {
    // set up a shear for stereo views
    if (this->LeftEye)
    {
      this->ProjectionTransform->Stereo(-this->EyeAngle / 2, this->Distance);
    }
    else
    {
      this->ProjectionTransform->Stereo(+this->EyeAngle / 2, this->Distance);
    }
  }

  if (this->ViewShear[0] != 0.0 || this->ViewShear[1] != 0.0)
  {
    this->ProjectionTransform->Shear(
      this->ViewShear[0], this->ViewShear[1], this->ViewShear[2] * this->Distance);
  }
}

//------------------------------------------------------------------------------
// Return the projection transform matrix. See ComputeProjectionTransform.
vtkMatrix4x4* vtkCamera::GetProjectionTransformMatrix(vtkRenderer* ren)
{
  double aspect[2];
  int lowerLeft[2];
  int usize, vsize;
  vtkMatrix4x4* matrix = vtkMatrix4x4::New();

  ren->GetTiledSizeAndOrigin(&usize, &vsize, lowerLeft, lowerLeft + 1);

  // some renderer subclasses may have more complicated computations for the
  // aspect ratio. So take that into account by computing the difference
  // between our simple aspect ratio and what the actual renderer is reporting.
  ren->ComputeAspect();
  ren->GetAspect(aspect);
  double aspect2[2];
  ren->vtkViewport::ComputeAspect();
  ren->vtkViewport::GetAspect(aspect2);
  double aspectModification = aspect[0] * aspect2[1] / (aspect[1] * aspect2[0]);

  if (usize && vsize)
  {
    matrix->DeepCopy(this->GetProjectionTransformMatrix(aspectModification * usize / vsize, -1, 1));
    matrix->Transpose();
  }

  return matrix;
}

//------------------------------------------------------------------------------
// Return the projection transform matrix. See ComputeProjectionTransform.
vtkMatrix4x4* vtkCamera::GetProjectionTransformMatrix(double aspect, double nearz, double farz)
{
  this->ComputeProjectionTransform(aspect, nearz, farz);

  // return the transform
  return this->ProjectionTransform->GetMatrix();
}

//------------------------------------------------------------------------------
// Return the projection transform object. See ComputeProjectionTransform.
vtkPerspectiveTransform* vtkCamera::GetProjectionTransformObject(
  double aspect, double nearz, double farz)
{
  this->ComputeProjectionTransform(aspect, nearz, farz);

  // return the transform
  return this->ProjectionTransform;
}

//------------------------------------------------------------------------------
// Return the projection transform matrix. See ComputeProjectionTransform.
vtkMatrix4x4* vtkCamera::GetCompositeProjectionTransformMatrix(
  double aspect, double nearz, double farz)
{
  // turn off stereo, the CompositeProjectionTransformMatrix is used for
  // picking, not for rendering.
  int stereo = this->Stereo;
  this->Stereo = 0;

  vtkMatrix4x4* projectionMatrix = this->GetProjectionTransformMatrix(aspect, nearz, farz);
  vtkMatrix4x4* viewMatrix = this->GetViewTransformMatrix();

  this->Transform->Identity();
  this->Transform->Concatenate(projectionMatrix);
  this->Transform->Concatenate(viewMatrix);

  this->Stereo = stereo;

  // return the transform
  return this->Transform->GetMatrix();
}

//------------------------------------------------------------------------------
// Return the attached light transform matrix.
vtkMatrix4x4* vtkCamera::GetCameraLightTransformMatrix()
{
  // return the transform
  return this->CameraLightTransform->GetMatrix();
}

//------------------------------------------------------------------------------
void vtkCamera::ComputeViewPlaneNormal()
{
  if (this->ViewShear[0] != 0.0 || this->ViewShear[1] != 0.0)
  {
    // set the VPN in camera coordinates
    this->ViewPlaneNormal[0] = this->ViewShear[0];
    this->ViewPlaneNormal[1] = this->ViewShear[1];
    this->ViewPlaneNormal[2] = 1.0;
    // transform the VPN to world coordinates using inverse of view transform
    this->ViewTransform->GetLinearInverse()->TransformNormal(
      this->ViewPlaneNormal, this->ViewPlaneNormal);
  }
  else
  {
    // VPN is -DOP
    this->ViewPlaneNormal[0] = -this->DirectionOfProjection[0];
    this->ViewPlaneNormal[1] = -this->DirectionOfProjection[1];
    this->ViewPlaneNormal[2] = -this->DirectionOfProjection[2];
  }
}

//------------------------------------------------------------------------------
// Return the 6 planes (Ax + By + Cz + D = 0) that bound
// the view frustum.
void vtkCamera::GetFrustumPlanes(double aspect, double planes[24])
{
  int i;
  double f, normals[6][4], matrix[4][4];

  // set up the normals
  for (i = 0; i < 6; i++)
  {
    normals[i][0] = 0.0;
    normals[i][1] = 0.0;
    normals[i][2] = 0.0;
    normals[i][3] = 1.0;
    // if i is even set to 1, if odd set to -1
    normals[i][i / 2] = 1 - (i % 2) * 2;
  }

  if (this->UseExplicitAspectRatio)
  {
    aspect = this->ExplicitAspectRatio;
  }

  // get the composite perspective matrix
  vtkMatrix4x4::DeepCopy(*matrix, this->GetCompositeProjectionTransformMatrix(aspect, -1, +1));

  // transpose the matrix for use with normals
  vtkMatrix4x4::Transpose(*matrix, *matrix);

  // transform the normals to world coordinates
  for (i = 0; i < 6; i++)
  {
    vtkMatrix4x4::MultiplyPoint(*matrix, normals[i], normals[i]);

    f = 1.0 /
      sqrt(normals[i][0] * normals[i][0] + normals[i][1] * normals[i][1] +
        normals[i][2] * normals[i][2]);

    planes[4 * i + 0] = normals[i][0] * f;
    planes[4 * i + 1] = normals[i][1] * f;
    planes[4 * i + 2] = normals[i][2] * f;
    planes[4 * i + 3] = normals[i][3] * f;
  }
}

void vtkCamera::UpdateIdealShiftScale(double aspect)
{
  double matrix[4][4];
  double imatrix[4][4];

  if (this->UseExplicitAspectRatio)
  {
    aspect = this->ExplicitAspectRatio;
  }

  // get the composite perspective matrix
  vtkMatrix4x4::DeepCopy(*matrix, this->GetCompositeProjectionTransformMatrix(aspect, -1, +1));
  vtkMatrix4x4::Invert(*matrix, *imatrix);

  double tmp[4];
  tmp[0] = 0;
  tmp[1] = 0;
  tmp[2] = -1;
  tmp[3] = 1;
  vtkMatrix4x4::MultiplyPoint(*imatrix, tmp, tmp);

  double shift[3];
  shift[0] = tmp[0] / tmp[3];
  shift[1] = tmp[1] / tmp[3];
  shift[2] = tmp[2] / tmp[3];

  tmp[0] = 1;
  tmp[1] = 1;
  tmp[2] = -1;
  tmp[3] = 1;
  vtkMatrix4x4::MultiplyPoint(*imatrix, tmp, tmp);

  tmp[0] /= tmp[3];
  tmp[1] /= tmp[3];
  tmp[2] /= tmp[3];

  double scale = sqrt(vtkMath::Distance2BetweenPoints(tmp, shift));

  // now snap
  if (fabs(log10(scale / this->NearPlaneScale)) > this->ShiftScaleThreshold)
  {
    this->NearPlaneScale = scale;
  }

  // our metric for shifting depends on scale
  double dist2 = vtkMath::Distance2BetweenPoints(this->NearPlaneShift, shift);
  if (dist2 && log10(sqrt(dist2) / this->NearPlaneScale) > this->ShiftScaleThreshold)
  {
    this->NearPlaneShift[0] = shift[0];
    this->NearPlaneShift[1] = shift[1];
    this->NearPlaneShift[2] = shift[2];
  }

  // now the focal point calcs
  tmp[0] = this->FocalPoint[0];
  tmp[1] = this->FocalPoint[1];
  tmp[2] = this->FocalPoint[2];
  tmp[3] = 1.0;
  vtkMatrix4x4::MultiplyPoint(*matrix, tmp, tmp);

  tmp[0] = 0.0;
  tmp[1] = 0.0;
  tmp[2] /= tmp[3];
  double fpdepth = tmp[2];
  tmp[3] = 1.0;
  vtkMatrix4x4::MultiplyPoint(*imatrix, tmp, tmp);

  shift[0] = tmp[0] / tmp[3];
  shift[1] = tmp[1] / tmp[3];
  shift[2] = tmp[2] / tmp[3];

  tmp[0] = 1;
  tmp[1] = 1;
  tmp[2] = fpdepth;
  tmp[3] = 1;
  vtkMatrix4x4::MultiplyPoint(*imatrix, tmp, tmp);

  tmp[0] /= tmp[3];
  tmp[1] /= tmp[3];
  tmp[2] /= tmp[3];

  scale = sqrt(vtkMath::Distance2BetweenPoints(tmp, shift));

  // now snap
  if (fabs(log10(scale / this->FocalPointScale)) > this->ShiftScaleThreshold)
  {
    this->FocalPointScale = scale;
  }

  // our metric for shifting depends on scale
  dist2 = vtkMath::Distance2BetweenPoints(this->FocalPointShift, shift);
  if (dist2 && log10(sqrt(dist2) / this->FocalPointScale) > this->ShiftScaleThreshold)
  {
    this->FocalPointShift[0] = shift[0];
    this->FocalPointShift[1] = shift[1];
    this->FocalPointShift[2] = shift[2];
  }
}

//------------------------------------------------------------------------------
vtkMTimeType vtkCamera::GetViewingRaysMTime()
{
  return this->ViewingRaysMTime.GetMTime();
}

//------------------------------------------------------------------------------
void vtkCamera::ViewingRaysModified()
{
  this->ViewingRaysMTime.Modified();
}

//------------------------------------------------------------------------------
// Description:
// Copy the properties of `source' into `this'.
// Copy pointers of matrices.
// \pre source_exists!=0
// \pre not_this: source!=this
void vtkCamera::ShallowCopy(vtkCamera* source)
{
  assert("pre: source_exists" && source != nullptr);
  assert("pre: not_this" && source != this);

  this->PartialCopy(source);

  // Shallow copy of matrices:
  if (this->UserTransform != nullptr)
  {
    this->UserTransform->Delete();
  }
  this->UserTransform = source->UserTransform;
  if (this->UserTransform != nullptr)
  {
    this->UserTransform->Register(this);
  }
  if (this->UserViewTransform != nullptr)
  {
    this->UserViewTransform->Delete();
  }
  this->UserViewTransform = source->UserViewTransform;
  if (this->UserViewTransform != nullptr)
  {
    this->UserViewTransform->Register(this);
  }

  if (this->ViewTransform != nullptr)
  {
    this->ViewTransform->Delete();
  }
  this->ViewTransform = source->ViewTransform;
  if (this->ViewTransform != nullptr)
  {
    this->ViewTransform->Register(this);
  }

  if (this->ProjectionTransform != nullptr)
  {
    this->ProjectionTransform->Delete();
  }
  this->ProjectionTransform = source->ProjectionTransform;
  if (this->ProjectionTransform != nullptr)
  {
    this->ProjectionTransform->Register(this);
  }

  if (this->Transform != nullptr)
  {
    this->Transform->Delete();
  }
  this->Transform = source->Transform;
  if (this->Transform != nullptr)
  {
    this->Transform->Register(this);
  }

  if (this->CameraLightTransform != nullptr)
  {
    this->CameraLightTransform->Delete();
  }
  this->CameraLightTransform = source->CameraLightTransform;
  if (this->CameraLightTransform != nullptr)
  {
    this->CameraLightTransform->Register(this);
  }

  if (this->EyeTransformMatrix != nullptr)
  {
    this->EyeTransformMatrix->Delete();
  }
  this->EyeTransformMatrix = source->EyeTransformMatrix;
  if (this->EyeTransformMatrix != nullptr)
  {
    this->EyeTransformMatrix->Register(this);
  }

  if (this->ProjectionPlaneOrientationMatrix != nullptr)
  {
    this->ProjectionPlaneOrientationMatrix->Delete();
  }
  this->ProjectionPlaneOrientationMatrix = source->ProjectionPlaneOrientationMatrix;
  if (this->ProjectionPlaneOrientationMatrix != nullptr)
  {
    this->ProjectionPlaneOrientationMatrix->Register(this);
  }

  if (this->ModelTransformMatrix != nullptr)
  {
    this->ModelTransformMatrix->Delete();
  }
  this->ModelTransformMatrix = source->ModelTransformMatrix;
  if (this->ModelTransformMatrix != nullptr)
  {
    this->ModelTransformMatrix->Register(this);
  }

  if (this->ModelViewTransform != nullptr)
  {
    this->ModelViewTransform->Delete();
  }
  this->ModelViewTransform = source->ModelViewTransform;
  if (this->ModelViewTransform != nullptr)
  {
    this->ModelViewTransform->Register(this);
  }
}

//------------------------------------------------------------------------------
// Description:
// Copy the properties of `source' into `this'.
// Copy the contents of the matrices.
// \pre source_exists!=0
// \pre not_this: source!=this
void vtkCamera::DeepCopy(vtkCamera* source)
{
  assert("pre: source_exists" && source != nullptr);
  assert("pre: not_this" && source != this);

  this->PartialCopy(source);

  // Deep copy the matrices:
  if (source->UserTransform == nullptr)
  {
    if (this->UserTransform != nullptr)
    {
      this->UserTransform->UnRegister(this);
      this->UserTransform = nullptr;
    }
  }
  else
  {
    if (this->UserTransform == nullptr)
    {
      this->UserTransform =
        static_cast<vtkHomogeneousTransform*>(source->UserTransform->MakeTransform());
    }
    this->UserTransform->DeepCopy(source->UserTransform);
  }

  if (source->UserViewTransform == nullptr)
  {
    if (this->UserViewTransform != nullptr)
    {
      this->UserViewTransform->UnRegister(this);
      this->UserViewTransform = nullptr;
    }
  }
  else
  {
    if (this->UserViewTransform == nullptr)
    {
      this->UserViewTransform =
        static_cast<vtkHomogeneousTransform*>(source->UserViewTransform->MakeTransform());
    }
    this->UserViewTransform->DeepCopy(source->UserViewTransform);
  }

  if (source->ViewTransform == nullptr)
  {
    if (this->ViewTransform != nullptr)
    {
      this->ViewTransform->UnRegister(this);
      this->ViewTransform = nullptr;
    }
  }
  else
  {
    if (this->ViewTransform == nullptr)
    {
      this->ViewTransform = static_cast<vtkTransform*>(source->ViewTransform->MakeTransform());
    }
    this->ViewTransform->DeepCopy(source->ViewTransform);
  }

  if (source->ProjectionTransform == nullptr)
  {
    if (this->ProjectionTransform != nullptr)
    {
      this->ProjectionTransform->UnRegister(this);
      this->ProjectionTransform = nullptr;
    }
  }
  else
  {
    if (this->ProjectionTransform == nullptr)
    {
      this->ProjectionTransform =
        static_cast<vtkPerspectiveTransform*>(source->ProjectionTransform->MakeTransform());
    }
    this->ProjectionTransform->DeepCopy(source->ProjectionTransform);
  }

  if (source->Transform == nullptr)
  {
    if (this->Transform != nullptr)
    {
      this->Transform->UnRegister(this);
      this->Transform = nullptr;
    }
  }
  else
  {
    if (this->Transform == nullptr)
    {
      this->Transform = static_cast<vtkPerspectiveTransform*>(source->Transform->MakeTransform());
    }
    this->Transform->DeepCopy(source->Transform);
  }

  if (source->CameraLightTransform == nullptr)
  {
    if (this->CameraLightTransform != nullptr)
    {
      this->CameraLightTransform->UnRegister(this);
      this->CameraLightTransform = nullptr;
    }
  }
  else
  {
    if (this->CameraLightTransform == nullptr)
    {
      this->CameraLightTransform =
        static_cast<vtkTransform*>(source->CameraLightTransform->MakeTransform());
    }
    this->CameraLightTransform->DeepCopy(source->CameraLightTransform);
  }

  if (source->ModelViewTransform == nullptr)
  {
    if (this->ModelViewTransform != nullptr)
    {
      this->ModelViewTransform->UnRegister(this);
      this->ModelViewTransform = nullptr;
    }
  }
  else
  {
    if (this->ModelViewTransform == nullptr)
    {
      this->ModelViewTransform =
        static_cast<vtkTransform*>(source->ModelViewTransform->MakeTransform());
    }
    this->ModelViewTransform->DeepCopy(source->ModelViewTransform);
  }

  if (source->ModelTransformMatrix == nullptr)
  {
    if (this->ModelTransformMatrix != nullptr)
    {
      this->ModelTransformMatrix->UnRegister(this);
      this->ModelTransformMatrix = nullptr;
    }
  }
  else
  {
    if (this->ModelTransformMatrix == nullptr)
    {
      this->ModelTransformMatrix =
        static_cast<vtkMatrix4x4*>(source->ModelTransformMatrix->NewInstance());
    }
    this->ModelTransformMatrix->DeepCopy(source->ModelTransformMatrix);
  }

  if (source->EyeTransformMatrix == nullptr)
  {
    if (this->EyeTransformMatrix != nullptr)
    {
      this->EyeTransformMatrix->UnRegister(this);
      this->EyeTransformMatrix = nullptr;
    }
  }
  else
  {
    if (this->EyeTransformMatrix == nullptr)
    {
      this->EyeTransformMatrix =
        static_cast<vtkMatrix4x4*>(source->EyeTransformMatrix->NewInstance());
    }
    this->EyeTransformMatrix->DeepCopy(source->EyeTransformMatrix);
  }

  if (source->ProjectionPlaneOrientationMatrix == nullptr)
  {
    if (this->ProjectionPlaneOrientationMatrix != nullptr)
    {
      this->ProjectionPlaneOrientationMatrix->UnRegister(this);
      this->ProjectionPlaneOrientationMatrix = nullptr;
    }
  }
  else
  {
    if (this->ProjectionPlaneOrientationMatrix == nullptr)
    {
      this->ProjectionPlaneOrientationMatrix =
        static_cast<vtkMatrix4x4*>(source->ProjectionPlaneOrientationMatrix->NewInstance());
    }
    this->ProjectionPlaneOrientationMatrix->DeepCopy(source->ProjectionPlaneOrientationMatrix);
  }
}

//------------------------------------------------------------------------------
// Description:
// Copy the ivars. Do nothing for the matrices.
// Called by ShallowCopy() and DeepCopy()
// \pre source_exists!=0
// \pre not_this: source!=this
void vtkCamera::PartialCopy(vtkCamera* source)
{
  assert("pre: source_exists" && source != nullptr);
  assert("pre: not_this" && source != this);

  int i;

  i = 0;
  while (i < 2)
  {
    this->WindowCenter[i] = source->WindowCenter[i];
    this->ObliqueAngles[i] = source->ObliqueAngles[i];
    this->ClippingRange[i] = source->ClippingRange[i];
    ++i;
  }
  i = 0;
  while (i < 3)
  {
    this->FocalPoint[i] = source->FocalPoint[i];
    this->Position[i] = source->Position[i];
    this->ViewUp[i] = source->ViewUp[i];
    this->DirectionOfProjection[i] = source->DirectionOfProjection[i];
    this->ViewPlaneNormal[i] = source->ViewPlaneNormal[i];
    this->ViewShear[i] = source->ViewShear[i];

    this->ScreenBottomLeft[i] = source->ScreenBottomLeft[i];
    this->ScreenBottomRight[i] = source->ScreenBottomRight[i];
    this->ScreenTopRight[i] = source->ScreenTopRight[i];
    this->ScreenCenter[i] = source->ScreenCenter[i];
    ++i;
  }

  this->ViewAngle = source->ViewAngle;
  this->EyeAngle = source->EyeAngle;
  this->ParallelProjection = source->ParallelProjection;
  this->ParallelScale = source->ParallelScale;
  this->Stereo = source->Stereo;
  this->LeftEye = source->LeftEye;
  this->Thickness = source->Thickness;
  this->Distance = source->Distance;
  this->UseHorizontalViewAngle = source->UseHorizontalViewAngle;
  this->UseOffAxisProjection = source->UseOffAxisProjection;
  this->OffAxisClippingAdjustment = source->OffAxisClippingAdjustment;

  this->FocalDisk = source->FocalDisk;
  this->FocalDistance = source->FocalDistance;
  this->EyeSeparation = source->EyeSeparation;

  this->ViewingRaysMTime = source->ViewingRaysMTime;
}

//------------------------------------------------------------------------------
void vtkCamera::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "ClippingRange: (" << this->ClippingRange[0] << ", " << this->ClippingRange[1]
     << ")\n";
  os << indent << "DirectionOfProjection: (" << this->DirectionOfProjection[0] << ", "
     << this->DirectionOfProjection[1] << ", " << this->DirectionOfProjection[2] << ")\n";
  os << indent << "Distance: " << this->Distance << "\n";
  os << indent << "EyeAngle: " << this->EyeAngle << "\n";
  os << indent << "FocalDisk: " << this->FocalDisk << "\n";
  os << indent << "FocalDistance: " << this->FocalDistance << "\n";
  os << indent << "FocalPoint: (" << this->FocalPoint[0] << ", " << this->FocalPoint[1] << ", "
     << this->FocalPoint[2] << ")\n";
  os << indent << "ViewShear: (" << this->ViewShear[0] << ", " << this->ViewShear[1] << ", "
     << this->ViewShear[2] << ")\n";
  os << indent << "ParallelProjection: " << (this->ParallelProjection ? "On\n" : "Off\n");
  os << indent << "ParallelScale: " << this->ParallelScale << "\n";
  os << indent << "Position: (" << this->Position[0] << ", " << this->Position[1] << ", "
     << this->Position[2] << ")\n";
  os << indent << "Stereo: " << (this->Stereo ? "On\n" : "Off\n");
  os << indent << "Left Eye: " << this->LeftEye << endl;
  os << indent << "Thickness: " << this->Thickness << "\n";
  os << indent << "ViewAngle: " << this->ViewAngle << "\n";
  os << indent << "UseHorizontalViewAngle: " << this->UseHorizontalViewAngle << "\n";
  os << indent << "UserTransform: ";
  if (this->UserTransform)
  {
    os << this->UserTransform << "\n";
  }
  else
  {
    os << "(none)\n";
  }
  if (this->UserViewTransform)
  {
    os << this->UserViewTransform << "\n";
  }
  else
  {
    os << "(none)\n";
  }
  os << indent << "FreezeFocalPoint: ";
  if (this->FreezeFocalPoint)
  {
    os << this->FreezeFocalPoint << "\n";
  }
  else
  {
    os << "(none)\n";
  }
  os << indent << "ViewPlaneNormal: (" << this->ViewPlaneNormal[0] << ", "
     << this->ViewPlaneNormal[1] << ", " << this->ViewPlaneNormal[2] << ")\n";
  os << indent << "ViewUp: (" << this->ViewUp[0] << ", " << this->ViewUp[1] << ", "
     << this->ViewUp[2] << ")\n";
  os << indent << "WindowCenter: (" << this->WindowCenter[0] << ", " << this->WindowCenter[1]
     << ")\n";

  os << indent << "UseOffAxisProjection: (" << this->UseOffAxisProjection << ")\n";

  os << indent << "ScreenBottomLeft: (" << this->ScreenBottomLeft[0] << ", "
     << this->ScreenBottomLeft[1] << ", " << this->ScreenBottomLeft[2] << ")\n";

  os << indent << "ScreenBottomRight: (" << this->ScreenBottomRight[0] << ", "
     << this->ScreenBottomRight[1] << ", " << this->ScreenBottomRight[2] << ")\n";

  os << indent << "ScreenTopRight: (" << this->ScreenTopRight[0] << ", " << this->ScreenTopRight[1]
     << ", " << this->ScreenTopRight[2] << ")\n";

  os << indent << "ScreenCenter: (" << this->ScreenCenter[0] << ", " << this->ScreenCenter[1]
     << ", " << this->ScreenCenter[2] << ")\n";

  os << indent << "OffAxisClippingAdjustment: (" << this->OffAxisClippingAdjustment << ")\n";
  os << indent << "EyeSeparation: (" << this->EyeSeparation << ")\n";

  os << indent << "ProjectionPlaneOrientationMatrix: (";
  if (this->ProjectionPlaneOrientationMatrix)
  {
    os << this->ProjectionPlaneOrientationMatrix << "\n";
    this->ProjectionPlaneOrientationMatrix->PrintSelf(os, indent.GetNextIndent());
  }
  else
  {
    os << "none";
  }
  os << indent << ")\n";

  os << indent << "EyeTransformMatrix: (" << this->EyeTransformMatrix << "\n";
  this->EyeTransformMatrix->PrintSelf(os, indent.GetNextIndent());
  os << indent << ")\n";

  os << indent << "ModelTransformMatrix: (" << this->ModelTransformMatrix << "\n";
  this->ModelTransformMatrix->PrintSelf(os, indent.GetNextIndent());
  os << indent << ")\n";

  os << indent << "ProjectionTransform: (" << this->ProjectionTransform << "\n";
  this->ProjectionTransform->PrintSelf(os, indent.GetNextIndent());
  os << indent << ")\n";
}

//------------------------------------------------------------------------------
void vtkCamera::SetEyePosition(double eyePosition[3])
{
  if (!eyePosition)
  {
    vtkErrorMacro(<< "ERROR: Invalid or nullptr eye position\n");
    return;
  }

  this->EyeTransformMatrix->SetElement(0, 3, eyePosition[0]);
  this->EyeTransformMatrix->SetElement(1, 3, eyePosition[1]);
  this->EyeTransformMatrix->SetElement(2, 3, eyePosition[2]);

  this->ComputeViewTransform();
  this->ComputeCameraLightTransform();

  this->Modified();
}

//------------------------------------------------------------------------------
void vtkCamera::GetEyePosition(double eyePosition[3])
{
  if (!eyePosition)
  {
    vtkErrorMacro(<< "ERROR: Invalid or nullptr eye position\n");
    return;
  }

  eyePosition[0] = this->EyeTransformMatrix->GetElement(0, 3);
  eyePosition[1] = this->EyeTransformMatrix->GetElement(1, 3);
  eyePosition[2] = this->EyeTransformMatrix->GetElement(2, 3);
}

//------------------------------------------------------------------------------
void vtkCamera::GetStereoEyePosition(double eyePosition[3])
{
  if (!eyePosition)
  {
    vtkErrorMacro(<< "ERROR: Invalid or nullptr eye position\n");
    return;
  }

  // Create an eye at the origin so it's easy to do the left/right shifting
  double E[4] = { 0.0, 0.0, 0.0, 1.0 };
  double shiftDistance = this->EyeSeparation / 2.0;

  if (this->LeftEye)
  {
    E[0] -= shiftDistance;
  }
  else
  {
    E[0] += shiftDistance;
  }

  // Now transform the "origin eye" to its real position and orientation
  this->EyeTransformMatrix->MultiplyPoint(E, E);
  eyePosition[0] = E[0];
  eyePosition[1] = E[1];
  eyePosition[2] = E[2];
}

//------------------------------------------------------------------------------
void vtkCamera::GetEyePlaneNormal(double normal[3])
{
  if (!normal)
  {
    vtkErrorMacro(<< "ERROR: Invalid or nullptr normal\n");
    return;
  }

  // Homogeneous normal.
  double localNormal[4];

  // In off-axis projection, the eye matrix is only used for position, and
  // the view direction is not considered, as the eye could be looking at
  // any arbitrary region of the screen.  So we approximate the eye plane
  // normal as the line originating at the center of the screen and pointing
  // at the eye position.
  localNormal[0] = this->EyeTransformMatrix->GetElement(0, 3) - this->ScreenCenter[0];
  localNormal[1] = this->EyeTransformMatrix->GetElement(1, 3) - this->ScreenCenter[1];
  localNormal[2] = this->EyeTransformMatrix->GetElement(2, 3) - this->ScreenCenter[2];
  localNormal[3] = 0.0;

  vtkMath::Normalize(localNormal);

  normal[0] = localNormal[0];
  normal[1] = localNormal[1];
  normal[2] = localNormal[2];
}

//------------------------------------------------------------------------------
vtkMatrix4x4* vtkCamera::GetModelViewTransformMatrix()
{
  this->ComputeModelViewMatrix();

  return this->ModelViewTransform->GetMatrix();
}

//------------------------------------------------------------------------------
vtkTransform* vtkCamera::GetModelViewTransformObject()
{
  this->ComputeModelViewMatrix();

  return this->ModelViewTransform;
}

//------------------------------------------------------------------------------
vtkMatrix4x4* vtkCamera::GetViewTransformMatrix()
{
  return this->GetModelViewTransformMatrix();
}

//------------------------------------------------------------------------------
vtkTransform* vtkCamera::GetViewTransformObject()
{
  return this->GetModelViewTransformObject();
}

//------------------------------------------------------------------------------
double* vtkCamera::GetOrientation()
{
  return this->ViewTransform->GetOrientation();
}

//------------------------------------------------------------------------------
double* vtkCamera::GetOrientationWXYZ()
{
  return this->ViewTransform->GetOrientationWXYZ();
}

//------------------------------------------------------------------------------
void vtkCamera::SetEyeTransformMatrix(const double elements[16])
{
  this->EyeTransformMatrix->Element[0][0] = elements[0];
  this->EyeTransformMatrix->Element[0][1] = elements[1];
  this->EyeTransformMatrix->Element[0][2] = elements[2];
  this->EyeTransformMatrix->Element[0][3] = elements[3];

  this->EyeTransformMatrix->Element[1][0] = elements[4];
  this->EyeTransformMatrix->Element[1][1] = elements[5];
  this->EyeTransformMatrix->Element[1][2] = elements[6];
  this->EyeTransformMatrix->Element[1][3] = elements[7];

  this->EyeTransformMatrix->Element[2][0] = elements[8];
  this->EyeTransformMatrix->Element[2][1] = elements[9];
  this->EyeTransformMatrix->Element[2][2] = elements[10];
  this->EyeTransformMatrix->Element[2][3] = elements[11];

  this->EyeTransformMatrix->Element[3][0] = elements[12];
  this->EyeTransformMatrix->Element[3][1] = elements[13];
  this->EyeTransformMatrix->Element[3][2] = elements[14];
  this->EyeTransformMatrix->Element[3][3] = elements[15];

  this->ComputeViewTransform();
  this->ComputeCameraLightTransform();

  this->Modified();
}

//------------------------------------------------------------------------------
void vtkCamera::SetModelTransformMatrix(const double elements[16])
{
  this->ModelTransformMatrix->Element[0][0] = elements[0];
  this->ModelTransformMatrix->Element[0][1] = elements[1];
  this->ModelTransformMatrix->Element[0][2] = elements[2];
  this->ModelTransformMatrix->Element[0][3] = elements[3];

  this->ModelTransformMatrix->Element[1][0] = elements[4];
  this->ModelTransformMatrix->Element[1][1] = elements[5];
  this->ModelTransformMatrix->Element[1][2] = elements[6];
  this->ModelTransformMatrix->Element[1][3] = elements[7];

  this->ModelTransformMatrix->Element[2][0] = elements[8];
  this->ModelTransformMatrix->Element[2][1] = elements[9];
  this->ModelTransformMatrix->Element[2][2] = elements[10];
  this->ModelTransformMatrix->Element[2][3] = elements[11];

  this->ModelTransformMatrix->Element[3][0] = elements[12];
  this->ModelTransformMatrix->Element[3][1] = elements[13];
  this->ModelTransformMatrix->Element[3][2] = elements[14];
  this->ModelTransformMatrix->Element[3][3] = elements[15];
  this->ModelTransformMatrix->Modified();
  this->Modified();
}
VTK_ABI_NAMESPACE_END
