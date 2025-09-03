/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2014-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
    Copyright (C) 2020-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "EigenMatrix.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class cmptType>
void Foam::EigenMatrix<cmptType>::tridiagonaliseSymmMatrix()
{
    for (label j = 0; j < n_; ++j)
    {
        EValsRe_[j] = EVecs_[n_ - 1][j];
    }

    // Householder reduction to tridiagonal form
    for (label i = n_ - 1; i > 0; --i)
    {
        // Scale to avoid under/overflow.
        cmptType scale = Zero;
        cmptType h = Zero;

        for (label k = 0; k < i; ++k)
        {
            scale = scale + mag(EValsRe_[k]);
        }

        if (scale == 0.0)
        {
            EValsIm_[i] = EValsRe_[i - 1];

            for (label j = 0; j < i; ++j)
            {
                EValsRe_[j] = EVecs_[i - 1][j];
                EVecs_[i][j] = Zero;
                EVecs_[j][i] = Zero;
            }
        }
        else
        {
            // Generate Householder vector
            for (label k = 0; k < i; ++k)
            {
                EValsRe_[k] /= scale;
                h += EValsRe_[k]*EValsRe_[k];
            }

            cmptType f = EValsRe_[i - 1];
            cmptType g = Foam::sqrt(h);

            if (f > 0)
            {
                g = -g;
            }

            EValsIm_[i] = scale*g;
            h -= f*g;
            EValsRe_[i - 1] = f - g;

            for (label j = 0; j < i; ++j)
            {
                EValsIm_[j] = Zero;
            }

            // Apply similarity transformation to remaining columns
            for (label j = 0; j < i; ++j)
            {
                f = EValsRe_[j];
                EVecs_[j][i] = f;
                g = EValsIm_[j] + EVecs_[j][j]*f;

                for (label k = j + 1; k <= i - 1; ++k)
                {
                    g += EVecs_[k][j]*EValsRe_[k];
                    EValsIm_[k] += EVecs_[k][j]*f;
                }

                EValsIm_[j] = g;
            }

            f = Zero;

            for (label j = 0; j < i; ++j)
            {
                EValsIm_[j] /= h;
                f += EValsIm_[j]*EValsRe_[j];
            }

            cmptType hh = f/(2.0*h);

            for (label j = 0; j < i; ++j)
            {
                EValsIm_[j] -= hh*EValsRe_[j];
            }

            for (label j = 0; j < i; ++j)
            {
                f = EValsRe_[j];
                g = EValsIm_[j];

                for (label k = j; k <= i - 1; ++k)
                {
                    EVecs_[k][j] -=
                            (f*EValsIm_[k] + g*EValsRe_[k]);
                }

                EValsRe_[j] = EVecs_[i - 1][j];
                EVecs_[i][j] = Zero;
            }
        }

        EValsRe_[i] = h;
    }

    // Accumulate transformations
    for (label i = 0; i < n_ - 1; ++i)
    {
        EVecs_[n_ - 1][i] = EVecs_[i][i];
        EVecs_[i][i] = cmptType(1);
        cmptType h = EValsRe_[i + 1];

        if (h != 0.0)
        {
            for (label k = 0; k <= i; ++k)
            {
                EValsRe_[k] = EVecs_[k][i + 1]/h;
            }

            for (label j = 0; j <= i; ++j)
            {
                cmptType g = Zero;

                for (label k = 0; k <= i; ++k)
                {
                    g += EVecs_[k][i + 1]*EVecs_[k][j];
                }

                for (label k = 0; k <= i; ++k)
                {
                    EVecs_[k][j] -= g*EValsRe_[k];
                }
            }
        }

        for (label k = 0; k <= i; ++k)
        {
            EVecs_[k][i + 1] = Zero;
        }
    }

    for (label j = 0; j < n_; ++j)
    {
        EValsRe_[j] = EVecs_[n_ - 1][j];
        EVecs_[n_ - 1][j] = Zero;
    }

    EVecs_[n_ - 1][n_ - 1] = cmptType(1);
    EValsIm_[0] = Zero;
}


template<class cmptType>
void Foam::EigenMatrix<cmptType>::symmTridiagonalQL()
{
    for (label i = 1; i < n_; ++i)
    {
        EValsIm_[i - 1] = EValsIm_[i];
    }

    EValsIm_[n_-1] = Zero;

    cmptType f = Zero;
    cmptType tst1 = Zero;
    cmptType eps = pow(2.0, -52.0);

    for (label l = 0; l < n_; l++)
    {
        // Find SMALL subdiagonal element
        tst1 = max(tst1, mag(EValsRe_[l]) + mag(EValsIm_[l]));
        label m = l;

        // Original while-loop from Java code
        while (m < n_)
        {
            if (mag(EValsIm_[m]) <= eps*tst1)
            {
                break;
            }

            ++m;
        }

        // If m == l, EValsRe_[l] is an eigenvalue, otherwise, iterate
        if (m > l)
        {
            do
            {
                // Compute implicit shift
                cmptType g = EValsRe_[l];
                cmptType p = (EValsRe_[l + 1] - g)/(2.0*EValsIm_[l]);
                cmptType r = Foam::hypot(p, cmptType(1));

                if (p < 0)
                {
                    r = -r;
                }

                EValsRe_[l] = EValsIm_[l]/(p + r);
                EValsRe_[l + 1] = EValsIm_[l]*(p + r);
                cmptType dl1 = EValsRe_[l + 1];
                cmptType h = g - EValsRe_[l];

                for (label i = l + 2; i < n_; ++i)
                {
                    EValsRe_[i] -= h;
                }

                f += h;

                // Implicit QL transformation.
                p = EValsRe_[m];
                cmptType c = cmptType(1);
                cmptType c2 = c;
                cmptType c3 = c;
                cmptType el1 = EValsIm_[l + 1];
                cmptType s = Zero;
                cmptType s2 = Zero;

                for (label i = m - 1; i >= l; --i)
                {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c*EValsIm_[i];
                    h = c*p;
                    r = std::hypot(p, EValsIm_[i]);
                    EValsIm_[i + 1] = s*r;
                    s = EValsIm_[i]/r;
                    c = p/r;
                    p = c*EValsRe_[i] - s*g;
                    EValsRe_[i + 1] = h + s*(c*g + s*EValsRe_[i]);

                    // Accumulate transformation
                    for (label k = 0; k < n_; ++k)
                    {
                        h = EVecs_[k][i + 1];
                        EVecs_[k][i + 1] = s*EVecs_[k][i] + c*h;
                        EVecs_[k][i] = c*EVecs_[k][i] - s*h;
                    }
                }

                p = -s*s2*c3*el1*EValsIm_[l]/dl1;
                EValsIm_[l] = s*p;
                EValsRe_[l] = c*p;
            }
            while (mag(EValsIm_[l]) > eps*tst1); // Convergence check
        }

        EValsRe_[l] = EValsRe_[l] + f;
        EValsIm_[l] = Zero;
    }

    // Sorting eigenvalues and corresponding vectors
    for (label i = 0; i < n_ - 1; ++i)
    {
        label k = i;
        cmptType p = EValsRe_[i];

        for (label j = i + 1; j < n_; ++j)
        {
            if (EValsRe_[j] < p)
            {
                k = j;
                p = EValsRe_[j];
            }
        }

        if (k != i)
        {
            EValsRe_[k] = EValsRe_[i];
            EValsRe_[i] = p;

            for (label j = 0; j < n_; ++j)
            {
                p = EVecs_[j][i];
                EVecs_[j][i] = EVecs_[j][k];
                EVecs_[j][k] = p;
            }
        }
    }
}


template<class cmptType>
void Foam::EigenMatrix<cmptType>::Hessenberg()
{
    DiagonalMatrix<cmptType> orth_(n_);

    label low = 0;
    label high = n_ - 1;

    for (label m = low + 1; m <= high - 1; ++m)
    {
        // Scale column
        cmptType scale = Zero;

        for (label i = m; i <= high; ++i)
        {
            scale = scale + mag(H_[i][m - 1]);
        }

        if (scale != 0.0)
        {
            // Compute Householder transformation
            cmptType h = Zero;

            for (label i = high; i >= m; --i)
            {
                orth_[i] = H_[i][m - 1]/scale;
                h += orth_[i]*orth_[i];
            }

            cmptType g = Foam::sqrt(h);

            if (orth_[m] > 0)
            {
                g = -g;
            }

            h -= orth_[m]*g;
            orth_[m] -= g;

            // Apply Householder similarity transformation
            // H = (I - u*u'/h)*H*(I - u*u')/h)
            for (label j = m; j < n_; ++j)
            {
                cmptType f = Zero;

                for (label i = high; i >= m; --i)
                {
                    f += orth_[i]*H_[i][j];
                }

                f /= h;

                for (label i = m; i <= high; ++i)
                {
                    H_[i][j] -= f*orth_[i];
                }
            }

            for (label i = 0; i <= high; ++i)
            {
                cmptType f = Zero;

                for (label j = high; j >= m; --j)
                {
                    f += orth_[j]*H_[i][j];
                }

                f /= h;

                for (label j = m; j <= high; ++j)
                {
                    H_[i][j] -= f*orth_[j];
                }
            }

            orth_[m] = scale*orth_[m];
            H_[m][m-1] = scale*g;
        }
    }

    // Accumulate transformations
    for (label i = 0; i < n_; ++i)
    {
        for (label j = 0; j < n_; ++j)
        {
            EVecs_[i][j] = (i == j ? 1.0 : 0.0);
        }
    }

    for (label m = high - 1; m >= low + 1; --m)
    {
        if (H_[m][m - 1] != 0.0)
        {
            for (label i = m + 1; i <= high; ++i)
            {
                orth_[i] = H_[i][m - 1];
            }

            for (label j = m; j <= high; ++j)
            {
                cmptType g = Zero;

                for (label i = m; i <= high; ++i)
                {
                    g += orth_[i]*EVecs_[i][j];
                }

                // Double division avoids possible underflow
                g = (g/orth_[m])/H_[m][m - 1];

                for (label i = m; i <= high; ++i)
                {
                    EVecs_[i][j] += g*orth_[i];
                }
            }
        }
    }
}


template<class cmptType>
void Foam::EigenMatrix<cmptType>::realSchur()
{
    // Initialize
    label nn = n_;
    label n = nn - 1;
    label low = 0;
    label high = nn - 1;
    cmptType eps = pow(2.0, -52.0);
    cmptType exshift = Zero;
    cmptType p = 0, q = 0, r = 0, s = 0, z = 0, t, w, x, y;

    // Store roots isolated by balance and compute matrix norm
    cmptType norm = Zero;

    for (label i = 0; i < nn; ++i)
    {
        if ((i < low) || (i > high))
        {
            EValsRe_[i] = H_[i][i];
            EValsIm_[i] = Zero;
        }

        for (label j = max(i - 1, 0); j < nn; ++j)
        {
            norm += mag(H_[i][j]);
        }
    }

    label iter = 0;

    while (n >= low)
    {
        // Look for single SMALL sub-diagonal element
        label l = n;

        while (l > low)
        {
            s = mag(H_[l - 1][l - 1]) + mag(H_[l][l]);

            if (s == 0.0)
            {
                s = norm;
            }

            if (mag(H_[l][l - 1]) < eps*s)
            {
                break;
            }

            l--;
        }

        // Check for convergence
        if (l == n)         // One root found
        {
            H_[n][n] = H_[n][n] + exshift;
            EValsRe_[n] = H_[n][n];
            EValsIm_[n] = Zero;
            --n;
            iter = 0;
        }
        else if (l == n - 1)  // Two roots found
        {
            w = H_[n][n - 1]*H_[n - 1][n];
            p = (H_[n - 1][n - 1] - H_[n][n])/2.0;
            q = p*p + w;
            z = Foam::sqrt(mag(q));
            H_[n][n] += exshift;
            H_[n - 1][n - 1] += exshift;
            x = H_[n][n];

            // Scalar pair
            if (q >= 0)
            {
                if (p >= 0)
                {
                    z = p + z;
                }
                else
                {
                    z = p - z;
                }

                EValsRe_[n - 1] = x + z;
                EValsRe_[n] = EValsRe_[n - 1];

                if (z != 0.0)
                {
                    EValsRe_[n] = x - w/z;
                }

                EValsIm_[n - 1] = Zero;
                EValsIm_[n] = Zero;
                x = H_[n][n - 1];
                s = mag(x) + mag(z);
                p = x/s;
                q = z/s;
                r = Foam::sqrt(p*p + q*q);
                p /= r;
                q /= r;

                // Row modification
                for (label j = n - 1; j < nn; ++j)
                {
                    z = H_[n - 1][j];
                    H_[n - 1][j] = q*z + p*H_[n][j];
                    H_[n][j] = q*H_[n][j] - p*z;
                }

                // Column modification
                for (label i = 0; i <= n; ++i)
                {
                    z = H_[i][n - 1];
                    H_[i][n - 1] = q*z + p*H_[i][n];
                    H_[i][n] = q*H_[i][n] - p*z;
                }

                // Accumulate transformations
                for (label i = low; i <= high; ++i)
                {
                    z = EVecs_[i][n - 1];
                    EVecs_[i][n - 1] = q*z + p*EVecs_[i][n];
                    EVecs_[i][n] = q*EVecs_[i][n] - p*z;
                }
            }
            else         // Complex pair of roots
            {
                EValsRe_[n - 1] = x + p;
                EValsRe_[n] = x + p;
                EValsIm_[n - 1] = z;
                EValsIm_[n] = -z;
            }

            n -= 2;
            iter = 0;
        }
        else            // No convergence yet
        {
            // Form shift
            x = H_[n][n];
            y = Zero;
            w = Zero;

            if (l < n)
            {
                y = H_[n - 1][n - 1];
                w = H_[n][n - 1]*H_[n - 1][n];
            }

            // Wilkinson's original ad-hoc shift
            if (iter == 10)
            {
                exshift += x;
                for (label i = low; i <= n; ++i)
                {
                    H_[i][i] -= x;
                }

                s = mag(H_[n][n - 1]) + mag(H_[n - 1][n - 2]);
                x = 0.75*s;
                y = 0.75*s;
                w = -0.4375*s*s;
            }

            // New ad hoc shift
            if (iter == 30)
            {
                s = (y - x)/2.0;
                s = s*s + w;

                if (s > 0)
                {
                    s = Foam::sqrt(s);

                    if (y < x)
                    {
                        s = -s;
                    }

                    s = x - w/((y - x)/2.0 + s);

                    for (label i = low; i <= n; ++i)
                    {
                        H_[i][i] -= s;
                    }

                    exshift += s;
                    x = 0.964;
                    y = 0.964;
                    w = 0.964;
                }
            }

            iter += 1;

            // Look for two consecutive SMALL sub-diagonal elements
            label m = n - 2;

            while (m >= l)
            {
                z = H_[m][m];
                r = x - z;
                s = y - z;
                p = (r*s - w)/H_[m + 1][m] + H_[m][m + 1];
                q = H_[m + 1][m + 1] - z - r - s;
                r = H_[m + 2][m + 1];
                s = mag(p) + mag(q) + mag(r);
                p /= s;
                q /= s;
                r /= s;

                if (m == l)
                {
                    break;
                }

                if
                (
                    mag(H_[m][m - 1])*(mag(q) + mag(r))
                            < eps*(mag(p)*(mag(H_[m - 1][m - 1])
                                    + mag(z) + mag(H_[m + 1][m + 1])))
                )
                {
                    break;
                }

                --m;
            }

            for (label i = m + 2; i <= n; ++i)
            {
                H_[i][i - 2] = Zero;

                if (i > m + 2)
                {
                    H_[i][i - 3] = Zero;
                }
            }

            // Double QR step involving rows l:n and columns m:n
            for (label k = m; k <= n - 1; ++k)
            {
                label notlast = (k != n - 1);

                if (k != m)
                {
                    p = H_[k][k - 1];
                    q = H_[k + 1][k - 1];
                    r = (notlast ? H_[k + 2][k - 1] : 0.0);
                    x = mag(p) + mag(q) + mag(r);

                    if (x != 0.0)
                    {
                        p /= x;
                        q /= x;
                        r /= x;
                    }
                }

                if (x == 0.0)
                {
                    break;
                }

                s = Foam::sqrt(p*p + q*q + r*r);

                if (p < 0)
                {
                    s = -s;
                }

                if (s != 0)
                {
                    if (k != m)
                    {
                        H_[k][k - 1] = -s*x;
                    }
                    else if (l != m)
                    {
                        H_[k][k - 1] = -H_[k][k - 1];
                    }

                    p = p + s;
                    x = p/s;
                    y = q/s;
                    z = r/s;
                    q /= p;
                    r /= p;

                    // Row modification
                    for (label j = k; j < nn; ++j)
                    {
                        p = H_[k][j] + q*H_[k + 1][j];

                        if (notlast)
                        {
                            p += r*H_[k + 2][j];
                            H_[k + 2][j] -= p*z;
                        }

                        H_[k][j] -= p*x;
                        H_[k + 1][j] -= p*y;
                    }

                    // Column modification
                    for (label i = 0; i <= min(n, k + 3); ++i)
                    {
                        p = x*H_[i][k] + y*H_[i][k + 1];

                        if (notlast)
                        {
                            p += z*H_[i][k + 2];
                            H_[i][k + 2] -= p*r;
                        }

                        H_[i][k] -= p;
                        H_[i][k + 1] -= p*q;
                    }

                    // Accumulate transformations
                    for (label i = low; i <= high; ++i)
                    {
                        p = x*EVecs_[i][k] + y*EVecs_[i][k + 1];

                        if (notlast)
                        {
                            p += z*EVecs_[i][k + 2];
                            EVecs_[i][k + 2] -= p*r;
                        }

                        EVecs_[i][k] -= p;
                        EVecs_[i][k + 1] -= p*q;
                    }
                }
            }
        }
    }

    // Backsubstitute to find vectors of upper triangular form
    if (norm == 0.0)
    {
        return;
    }

    for (n = nn-1; n >= 0; --n)
    {
        p = EValsRe_[n];
        q = EValsIm_[n];

        // scalar vector
        if (q == 0)
        {
            label l = n;
            H_[n][n] = cmptType(1);

            for (label i = n-1; i >= 0; --i)
            {
                w = H_[i][i] - p;
                r = Zero;

                for (label j = l; j <= n; ++j)
                {
                    r += H_[i][j]*H_[j][n];
                }

                if (EValsIm_[i] < 0.0)
                {
                    z = w;
                    s = r;
                }
                else
                {
                    l = i;

                    if (EValsIm_[i] == 0.0)
                    {
                        if (w != 0.0)
                        {
                            H_[i][n] = -r/w;
                        }
                        else
                        {
                            H_[i][n] = -r/(eps*norm);
                        }
                    }
                    else // Solve real equations
                    {
                        x = H_[i][i + 1];
                        y = H_[i + 1][i];
                        q = (EValsRe_[i] - p)*(EValsRe_[i] - p)
                                + EValsIm_[i]*EValsIm_[i];

                        t = (x*s - z*r)/q;
                        H_[i][n] = t;

                        if (mag(x) > mag(z))
                        {
                            H_[i + 1][n] = (-r - w*t)/x;
                        }
                        else
                        {
                            H_[i + 1][n] = (-s - y*t)/z;
                        }
                    }

                    // Overflow control
                    t = mag(H_[i][n]);

                    if ((eps*t)*t > 1)
                    {
                        for (label j = i; j <= n; ++j)
                        {
                            H_[j][n] /= t;
                        }
                    }
                }
            }
        }
        else if (q < 0) // Complex vector
        {
            label l = n - 1;

            // Last vector component imaginary so matrix is triangular
            if (mag(H_[n][n - 1]) > mag(H_[n - 1][n]))
            {
                H_[n - 1][n - 1] = q/H_[n][n - 1];
                H_[n - 1][n] = -(H_[n][n] - p)/H_[n][n - 1];
            }
            else
            {
                complex cDiv = complex(0.0, -H_[n - 1][n])
                        /complex(H_[n - 1][n-1]-p, q);

                H_[n - 1][n - 1] = cDiv.real();
                H_[n - 1][n] = cDiv.imag();
            }

            H_[n][n - 1] = Zero;
            H_[n][n] = cmptType(1);

            for (label i = n - 2; i >= 0; --i)
            {
                cmptType ra, sa, vr, vi;
                ra = Zero;
                sa = Zero;

                for (label j = l; j <= n; ++j)
                {
                    ra += H_[i][j]*H_[j][n-1];
                    sa += H_[i][j]*H_[j][n];
                }

                w = H_[i][i] - p;

                if (EValsIm_[i] < 0.0)
                {
                    z = w;
                    r = ra;
                    s = sa;
                }
                else
                {
                    l = i;

                    if (EValsIm_[i] == 0)
                    {
                        complex cDiv = complex(-ra, -sa)/complex(w, q);
                        H_[i][n - 1] = cDiv.real();
                        H_[i][n] = cDiv.imag();
                    }
                    else
                    {
                        // Solve complex equations
                        x = H_[i][i + 1];
                        y = H_[i + 1][i];
                        vr = (EValsRe_[i] - p)*(EValsRe_[i] - p)
                                + EValsIm_[i]*EValsIm_[i] - q*q;

                        vi = 2.0*(EValsRe_[i] - p)*q;

                        if ((vr == 0.0) && (vi == 0.0))
                        {
                            vr = eps*norm*(mag(w) + mag(q) + mag(x) + mag(y)
                                    + mag(z));
                        }

                        complex cDiv =
                                complex(x*r - z*ra + q*sa, x*s - z*sa - q*ra)
                                    /complex(vr, vi);

                        H_[i][n - 1] = cDiv.real();
                        H_[i][n] = cDiv.imag();

                        if (mag(x) > (mag(z) + mag(q)))
                        {
                            H_[i + 1][n - 1] = (-ra - w*H_[i][n - 1]
                                    + q*H_[i][n])/x;

                            H_[i + 1][n] = (-sa - w*H_[i][n]
                                    - q*H_[i][n - 1])/x;
                        }
                        else
                        {
                            complex cDiv = complex(-r - y*H_[i][n - 1], -s
                                    - y*H_[i][n])/complex(z, q);

                            H_[i + 1][n - 1] = cDiv.real();
                            H_[i + 1][n] = cDiv.imag();
                        }
                    }

                    // Overflow control
                    t = max(mag(H_[i][n - 1]), mag(H_[i][n]));

                    if ((eps*t)*t > 1)
                    {
                        for (label j = i; j <= n; ++j)
                        {
                            H_[j][n - 1] /= t;
                            H_[j][n] /= t;
                        }
                    }
                }
            }
        }
    }

    // Vectors of isolated roots
    for (label i = 0; i < nn; ++i)
    {
        if (i < low || i > high)
        {
            for (label j = i; j < nn; ++j)
            {
                EVecs_[i][j] = H_[i][j];
            }
        }
    }

    // Back transformation to get eigenvectors of original matrix
    for (label j = nn - 1; j >= low; --j)
    {
        for (label i = low; i <= high; ++i)
        {
            z = Zero;

            for (label k = low; k <= min(j, high); ++k)
            {
                z = z + EVecs_[i][k]*H_[k][j];
            }

            EVecs_[i][j] = z;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class cmptType>
Foam::EigenMatrix<cmptType>::EigenMatrix(const SquareMatrix<cmptType>& A)
:
    n_(A.n()),
    EValsRe_(n_, Zero),
    EValsIm_(n_, Zero),
    EVecs_(n_, Zero),
    H_()
{
    if (n_ <= 0)
    {
        FatalErrorInFunction
            << "Input matrix has zero size."
            << abort(FatalError);
    }

    if (A.symmetric())
    {
        EVecs_ = A;
        tridiagonaliseSymmMatrix();
        symmTridiagonalQL();
    }
    else
    {
        H_ = A;
        Hessenberg();
        realSchur();
    }
}


template<class cmptType>
Foam::EigenMatrix<cmptType>::EigenMatrix
(
    const SquareMatrix<cmptType>& A,
    bool symmetric
)
:
    n_(A.n()),
    EValsRe_(n_, Zero),
    EValsIm_(n_, Zero),
    EVecs_(n_, Zero),
    H_()
{
    if (n_ <= 0)
    {
        FatalErrorInFunction
            << "Input matrix has zero size."
            << abort(FatalError);
    }

    if (symmetric)
    {
        EVecs_ = A;
        tridiagonaliseSymmMatrix();
        symmTridiagonalQL();
    }
    else
    {
        H_ = A;
        Hessenberg();
        realSchur();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class cmptType>
const Foam::SquareMatrix<Foam::complex>
Foam::EigenMatrix<cmptType>::complexEVecs() const
{
    SquareMatrix<complex> EVecs(EVecs_.n());

    std::transform
    (
        EVecs_.cbegin(),
        EVecs_.cend(),
        EVecs.begin(),
        [&](const cmptType& x){ return complex(x); }
    );

    for (label i = 0; i < EValsIm_.size(); ++i)
    {
        if (mag(EValsIm_[i]) > VSMALL)
        {
            for (label j = 0; j < EVecs.m(); ++j)
            {
                EVecs(j, i) = complex(EVecs_(j, i), EVecs_(j, i+1));
                EVecs(j, i+1) = EVecs(j, i).conjugate();
            }
            ++i;
        }
    }

    return EVecs;
}


// ************************************************************************* //
