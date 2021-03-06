/*
 * Copyright (c) 2008-2011 Zhang Ming (M. Zhang), zmjerry@163.com
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 2 or any later version.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details. A copy of the GNU General Public License is available at:
 * http://www.fsf.org/licensing/licenses
 */


/*****************************************************************************
 *                             correlation-impl.h
 *
 * Implementation for linear correlation.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/
#include "correlation.h"

namespace splab {

/**
 * Auto-correlation by defination in time domain.
 */
    template<typename Type>
    inline Vector<Type> corr( const Vector<Type> &xn, const string &opt )
    {
        Vector<Type> rn = conv( xn, reverse(xn) );

        biasedProcessing( rn, opt );

        return rn;
    }


/**
 * Cross-correlation by defination in time domain.
 */
    template<typename Type>
    inline Vector<Type> corr( const Vector<Type> &xn, const Vector<Type> &yn,
                              const string &opt )
    {
        int N = xn.size(),
                d = N - yn.size();
        Vector<Type> rn;

        if( d > 0 )
            rn = conv( xn, reverse(wextend(yn,d,"right","zpd")) );
        else if( d < 0 )
        {
            N -= d;
            rn = conv( wextend(xn,-d,"right","zpd"), reverse(yn) );
        }
        else
            rn = conv( xn, reverse(yn) );

        biasedProcessing( rn, opt );

        return rn;
    }


/**
 * Fast auto-correlation by using FFT.
 */
    template<typename Type>
    inline Vector<Type> fastCorr( const Vector<Type> &xn, const string &opt )
    {
        Vector<Type> rn = fastConv( xn, reverse(xn) );

        biasedProcessing( rn, opt );

        return rn;
    }


/**
 * Fast cross-correlation by using FFT.
 */
    template<typename Type>
    inline Vector<Type> fastCorr( const Vector<Type> &xn, const Vector<Type> &yn,
                                  const string &opt )
    {
        int N = xn.size(),
                d = N - yn.size();
        Vector<Type> rn;

        if( d > 0 )
            rn = fastConv( xn, reverse(wextend(yn,d,"right","zpd")) );
        else if( d < 0 )
        {
            N -= d;
            rn = fastConv( wextend(xn,-d,"right","zpd"), reverse(yn) );
        }
        else
            rn = fastConv( xn, reverse(yn) );

        biasedProcessing( rn, opt );

        return rn;
    }


/**
 * Biase processing for correlation.
 */
    template<typename Type>
    static void biasedProcessing( Vector<Type> &rn, const string &opt )
    {
        int N = (rn.size()+1) / 2;

        if( opt == "biased" )
            rn /= Type(N);
        else if( opt == "unbiased" )
        {
            int mid = N-1;
            rn[mid] /= N;
            for( int i=1; i<N; ++i )
            {
                rn[mid+i] /= (N-i);
                rn[mid-i] /= (N-i);
            }
        }
    }
}