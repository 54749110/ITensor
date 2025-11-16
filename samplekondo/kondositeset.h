//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_SITESET_MORE_H
#define __ITENSOR_SITESET_MORE_H
#include "itensor/itensor.h"
#include "itensor/util/str.h"

namespace itensor {

//
// Classes derived from SiteSet 
// represent the Hilbert space of a 
// system as a set of Site indices.
//
// The convention for operators is
// that they are 2-index ITensors
// with the Site Index pointing
// In and the Site' Index pointing
// Out. This is so we can compute expectation
// values by doing dag(prime(A,Site)) * Op * A.
// (assuming the tensor A is an ortho center 
// of our MPS)
//

template<typename ASiteType, typename BSiteType, typename CSiteType, typename DSiteType, typename ESiteType>
class Tensor5SiteSet : public SiteSet
    {
    public:

    Tensor5SiteSet() { }

    Tensor5SiteSet(int N, 
                 Args const& args = Args::global())
        {
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
            {
            if(j%5 == 1) sites.set(j,ASiteType({args,"SiteNumber=",j}));
            else if(j%5 == 2) sites.set(j,BSiteType({args,"SiteNumber=",j}));
            else if(j%5 == 3) sites.set(j,CSiteType({args,"SiteNumber=",j}));
            else if(j%5 == 4) sites.set(j,DSiteType({args,"SiteNumber=",j}));
            else         sites.set(j,ESiteType({args,"SiteNumber=",j}));
            }
        SiteSet::init(std::move(sites));
        }

    Tensor5SiteSet(IndexSet const& is)
        {
        int N = is.length();
        auto sites = SiteStore(N);
        for(auto j : range1(N))
            {
            if(j%5 == 1) sites.set(j,ASiteType(is(j)));
            else if(j%5 == 2) sites.set(j,BSiteType(is(j)));
            else if(j%5 == 3) sites.set(j,CSiteType(is(j)));
            else if(j%5 == 4) sites.set(j,DSiteType(is(j)));
            else         sites.set(j,ESiteType(is(j)));
            }
        SiteSet::init(std::move(sites));
        }

    void
    read(std::istream& s)
        {
        int N = itensor::read<int>(s);
        if(N > 0)
            {
            auto store = SiteStore(N);
            for(int j = 1; j <= N; ++j) 
                {
                auto I = Index{};
                I.read(s);
                if(j%5==1) store.set(j,ASiteType(I));
                else if(j%5==2) store.set(j,BSiteType(I));
                else if(j%5==3) store.set(j,CSiteType(I));
                else if(j%5==4) store.set(j,DSiteType(I));
                else       store.set(j,ESiteType(I));
                }
            init(std::move(store));
            }
        }
    };

template<typename ASiteType, typename BSiteType, typename CSiteType, typename DSiteType>
class Tensor4SiteSet : public SiteSet
    {
    public:

    Tensor4SiteSet() { }

    Tensor4SiteSet(int N, 
                 Args const& args = Args::global())
        {
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
            {
            if(j%4 == 1) sites.set(j,ASiteType({args,"SiteNumber=",j}));
            else if(j%4 == 2) sites.set(j,BSiteType({args,"SiteNumber=",j}));
            else if(j%4 == 3) sites.set(j,CSiteType({args,"SiteNumber=",j}));
            else         sites.set(j,DSiteType({args,"SiteNumber=",j}));
            }
        SiteSet::init(std::move(sites));
        }

    Tensor4SiteSet(IndexSet const& is)
        {
        int N = is.length();
        auto sites = SiteStore(N);
        for(auto j : range1(N))
            {
            if(j%4 == 1) sites.set(j,ASiteType(is(j)));
            else if(j%4 == 2) sites.set(j,BSiteType(is(j)));
            else if(j%4 == 3) sites.set(j,CSiteType(is(j)));
            else         sites.set(j,DSiteType(is(j)));
            }
        SiteSet::init(std::move(sites));
        }

    void
    read(std::istream& s)
        {
        int N = itensor::read<int>(s);
        if(N > 0)
            {
            auto store = SiteStore(N);
            for(int j = 1; j <= N; ++j) 
                {
                auto I = Index{};
                I.read(s);
                if(j%4==1) store.set(j,ASiteType(I));
                else if(j%4==2) store.set(j,BSiteType(I));
                else if(j%4==3) store.set(j,CSiteType(I));
                else       store.set(j,DSiteType(I));
                }
            init(std::move(store));
            }
        }
    };

template<typename ASiteType, typename BSiteType, typename CSiteType>
class Kondo2ChannelSiteSet : public SiteSet
    {
    public:

    Kondo2ChannelSiteSet() { }

    Kondo2ChannelSiteSet(int N, 
                 Args const& args = Args::global())
        {
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
            {
            if(j%3 == 1) sites.set(j,ASiteType({args,"SiteNumber=",j}));
            else if(j%3 == 2) sites.set(j,BSiteType({args,"SiteNumber=",j}));
            else         sites.set(j,CSiteType({args,"SiteNumber=",j}));
            }
        SiteSet::init(std::move(sites));
        }

    Kondo2ChannelSiteSet(IndexSet const& is)
        {
        int N = is.length();
        auto sites = SiteStore(N);
        for(auto j : range1(N))
            {
            if(j%3 == 1) sites.set(j,ASiteType(is(j)));
            else if(j%3 == 2) sites.set(j,BSiteType(is(j)));
            else         sites.set(j,CSiteType(is(j)));
            }
        SiteSet::init(std::move(sites));
        }

    void
    read(std::istream& s)
        {
        int N = itensor::read<int>(s);
        if(N > 0)
            {
            auto store = SiteStore(N);
            for(int j = 1; j <= N; ++j) 
                {
                auto I = Index{};
                I.read(s);
                if(j%3==1) store.set(j,ASiteType(I));
                else if(j%3==2) store.set(j,BSiteType(I));
                else       store.set(j,CSiteType(I));
                }
            init(std::move(store));
            }
        }
    };

template<typename ASiteType, typename BSiteType, typename CSiteType>
class Kondo2ChannelImpuritySiteSet : public SiteSet
    {
    public:

    Kondo2ChannelImpuritySiteSet() { }

    Kondo2ChannelImpuritySiteSet(int N, 
                 Args const& args = Args::global())
        {
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
            {
            if(j%2 == 1 && j<N) sites.set(j,ASiteType({args,"SiteNumber=",j}));
            else if(j%2 == 0 && j<N) sites.set(j,BSiteType({args,"SiteNumber=",j}));
            else         sites.set(j,CSiteType({args,"SiteNumber=",j}));
            }
        SiteSet::init(std::move(sites));
        }

    Kondo2ChannelImpuritySiteSet(IndexSet const& is)
        {
        int N = is.length();
        auto sites = SiteStore(N);
        for(auto j : range1(N))
            {
            if(j%2 == 1 && j<N) sites.set(j,ASiteType(is(j)));
            else if(j%2 == 0 && j<N) sites.set(j,BSiteType(is(j)));
            else         sites.set(j,CSiteType(is(j)));
            }
        SiteSet::init(std::move(sites));
        }

    void
    read(std::istream& s)
        {
        int N = itensor::read<int>(s);
        if(N > 0)
            {
            auto store = SiteStore(N);
            for(int j = 1; j <= N; ++j) 
                {
                auto I = Index{};
                I.read(s);
                if(j%2==1 && j<N) store.set(j,ASiteType(I));
                else if(j%2==0 && j<N) store.set(j,BSiteType(I));
                else       store.set(j,CSiteType(I));
                }
            init(std::move(store));
            }
        }
    };

template<typename ASiteType, typename BSiteType>
class Kondo2ChannelSiteSetv2 : public SiteSet
    {
    public:

    Kondo2ChannelSiteSetv2() { }

    Kondo2ChannelSiteSetv2(int N, 
                 Args const& args = Args::global())
        {
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
            {
            if(j > 2*N/3 ) sites.set(j,BSiteType({args,"SiteNumber=",j}));
            else         sites.set(j,ASiteType({args,"SiteNumber=",j}));
            }
        SiteSet::init(std::move(sites));
        }

    Kondo2ChannelSiteSetv2(IndexSet const& is)
        {
        int N = is.length();
        auto sites = SiteStore(N);
        for(auto j : range1(N))
            {
            if(j > 2*N/3 ) sites.set(j,BSiteType(is(j)));
            else         sites.set(j,ASiteType(is(j)));
            }
        SiteSet::init(std::move(sites));
        }

    void
    read(std::istream& s)
        {
        int N = itensor::read<int>(s);
        if(N > 0)
            {
            auto store = SiteStore(N);
            for(int j = 1; j <= N; ++j) 
                {
                auto I = Index{};
                I.read(s);
                if(j>2*N/3) store.set(j,BSiteType(I));
                else        store.set(j,ASiteType(I));
                }
            init(std::move(store));
            }
        }
    };

template<typename ASiteType, typename BSiteType>
class Kondo2ChannelSiteSetv3 : public SiteSet
    {
    public:

    Kondo2ChannelSiteSetv3() { }

    Kondo2ChannelSiteSetv3(int N, 
                 Args const& args = Args::global())
        {
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
            {
            if(j%3 == 2) sites.set(j,BSiteType({args,"SiteNumber=",j}));
            else         sites.set(j,ASiteType({args,"SiteNumber=",j}));
            }
        SiteSet::init(std::move(sites));
        }

    Kondo2ChannelSiteSetv3(IndexSet const& is)
        {
        int N = is.length();
        auto sites = SiteStore(N);
        for(auto j : range1(N))
            {
            if(j%3 == 2) sites.set(j,BSiteType(is(j)));
            else         sites.set(j,ASiteType(is(j)));
            }
        SiteSet::init(std::move(sites));
        }

    void
    read(std::istream& s)
        {
        int N = itensor::read<int>(s);
        if(N > 0)
            {
            auto store = SiteStore(N);
            for(int j = 1; j <= N; ++j) 
                {
                auto I = Index{};
                I.read(s);
                if(j%3 == 2) store.set(j,BSiteType(I));
                else         store.set(j,ASiteType(I));
                }
            init(std::move(store));
            }
        }
    };

template<typename ASiteType, typename BSiteType>
class Kondo2ChannelImpuritySiteSetv2 : public SiteSet
    {
    public:

    Kondo2ChannelImpuritySiteSetv2() { }

    Kondo2ChannelImpuritySiteSetv2(int N, 
                 Args const& args = Args::global())
        {
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
            {
            if(j == N) sites.set(j,BSiteType({args,"SiteNumber=",j}));
            else       sites.set(j,ASiteType({args,"SiteNumber=",j}));
            }
        SiteSet::init(std::move(sites));
        }

    Kondo2ChannelImpuritySiteSetv2(IndexSet const& is)
        {
        int N = is.length();
        auto sites = SiteStore(N);
        for(auto j : range1(N))
            {
            if(j == N) sites.set(j,BSiteType(is(j)));
            else       sites.set(j,ASiteType(is(j)));
            }
        SiteSet::init(std::move(sites));
        }

    void
    read(std::istream& s)
        {
        int N = itensor::read<int>(s);
        if(N > 0)
            {
            auto store = SiteStore(N);
            for(int j = 1; j <= N; ++j) 
                {
                auto I = Index{};
                I.read(s);
                if(j == N) store.set(j,BSiteType(I));
                else       store.set(j,ASiteType(I));
                }
            init(std::move(store));
            }
        }
    };



template<typename ASiteType, typename BSiteType>
class twochannelkondo : public SiteSet
    {
    public:

    twochannelkondo() { }

    twochannelkondo(int N, 
                 Args const& args = Args::global())
        {
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
            {
            if((j == 1 || j==N )) sites.set(j,ASiteType({args,"SiteNumber=",j}));
            else                 sites.set(j,BSiteType({args,"SiteNumber=",j}));
            }
        SiteSet::init(std::move(sites));
        }

    twochannelkondo(IndexSet const& is)
        {
        int N = is.length();
        auto sites = SiteStore(N);
        for(auto j : range1(N))
            {
            if((j == 1 || j==N)) sites.set(j,ASiteType(is(j)));
            else         sites.set(j,BSiteType(is(j)));
            }
        SiteSet::init(std::move(sites));
        }

    void
    read(std::istream& s)
        {
        int N = itensor::read<int>(s);
        if(N > 0)
            {
            auto store = SiteStore(N);
            for(int j = 1; j <= N; ++j) 
                {
                auto I = Index{};
                I.read(s);
                if((j == 1 || j==N)) store.set(j,ASiteType(I));
                else       store.set(j,BSiteType(I));
                }
            init(std::move(store));
            }
        }
    };

} //namespace itensor

#endif
