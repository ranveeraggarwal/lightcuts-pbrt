
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// core/integrator.cpp*
#include "stdafx.h"
#include "integrator.h"
#include "scene.h"
#include "intersection.h"
#include "montecarlo.h"
#include "lighttree.h"
#include <vector>
#include <iostream>
#define THRESHOLD 0.02
using namespace std;
// Integrator Method Definitions
Integrator::~Integrator() {
}
struct cutvertex_struct {
    LightNode * l;
    Spectrum error;
    Spectrum illumination;
};
typedef struct cutvertex_struct cutvertex;

bool operator <(const cutvertex& a, const cutvertex& b) { // Compare cut vertices on the basis of error
    float cola[3];
    float colb[3];
    a.error.ToRGB(cola);
    b.error.ToRGB(colb);
    if(colb[0] < 0 || colb[1] < 0 || colb[2] < 0) // Negative value means infinite error
        return true;
    if(cola[0] + cola[1] + cola[2] < colb[0] + colb[1] + colb[2])
        return true;
    return false;
}

// Set zero error for leaf. Call ErrorBound for others.
#define SetError(v) \
    do { \
    if(v.l->isLeaf) \
        v.error = zero; \
    else { \
        v.error = ErrorBound(p, n, v.l->minBoundBox, v.l->maxBoundBox, v.l->mainLight->lightPos, bsdf->f(wo, wi, BxDFType(BSDF_ALL & ~BSDF_SPECULAR)), v.l->Intensity);\
    } \
    } while(0)

// Use PBRT's illumination computation
#define Illuminate(v) \
    do { \
        PointLight *light = v.l->mainLight; \
        Spectrum temp = light->Intensity; \
        light->Intensity = v.l->Intensity; \
        int nSamples = 1; \
        Spectrum Ld(0.); \
        for (int j = 0; j < nSamples; ++j) { \
            LightSample lightSample; \
            BSDFSample bsdfSample; \
            lightSample = LightSample(rng); \
            bsdfSample = BSDFSample(rng); \
            Ld += EstimateDirect(scene, renderer, arena, light, p, n, wo, \
                rayEpsilon, time, bsdf, rng, lightSample, bsdfSample, \
                BxDFType(BSDF_ALL & ~BSDF_SPECULAR)); \
        } \
        v.illumination =  Ld / nSamples; \
        light->Intensity = temp; \
    } while(0)

// Function for the upper bound on error
Spectrum ErrorBound(Point p, Normal n, Point bb_low, Point bb_high, Point l, Spectrum brdf, Spectrum intensities) {
    double err = 1;
    bool inf = false;
    int i;
    double dist = 0;
    // Calculate minimum distance from bounding box
    for(i = 0; i < 3; i++) {
        if(p[i] < bb_low[i])
            dist += (p[i] - bb_low[i])*(p[i] - bb_low[i]);
        else if(p[i] > bb_high[i])
            dist += (p[i] - bb_high[i])*(p[i] - bb_high[i]);
        else {
            inf = true;
            break;
        }
    }
    // Infinite error bound if p is inside the bounding box
    if(inf) {
        return -1 * intensities;
    }
    // Bound of the geometric term
    err *= 1.0 / dist;
    n /= n.Length();
    // To bound the cosine term, align normal to the z-axis 
    // and take maximum z-coordinate and minimum absolute values of x and y coordinates.
    float M[3][3]; //Rotation matrix from normal to z-axis
    float v[3][3];
    float vsq[3][3];
    int j;
    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            if(i == j)
                M[i][j] = 1;
            else
                M[i][j] = 0;
            v[i][j] = 0;
            vsq[i][j] = 0;
        }
    }
    float vx = n[1];
    float vy = -(n[0]);
    float s = sqrtf(vx*vx + vy*vy);
    float c = n[2];
    Point bb_low2, bb_high2;
    for(i = 0; i < 3; i++) {
        bb_low2[i] = bb_low[i] - p[i];
        bb_high2[i] = bb_high[i] - p[i];
    }
    int k;
    if(fabsf(s) < 0.00001) {
        if(c < 0) {
            M[2][2] = -1;
        }
    }
    else {
        v[0][2] = vy;
        M[0][2] += vy;
        v[1][2] = -vx;
        M[1][2] += -vx;
        v[2][0] = -vy;
        M[2][0] += -vy;
        v[2][1] = vx;
        M[2][1] += vx;
        float temp = (1 - c) / (s * s);
        for(i = 0; i < 3; i++) {
            for(j = 0; j < 3; j++) {
                for(k = 0; k < 3; k++) {
                    vsq[i][j] += v[i][k] * v[k][j];
                }
                M[i][j] += vsq[i][j] * temp;
            }
        }
    }
    // Rotation Matrix computed
    Point bb_low3, bb_high3;
    // Transform and align normal to the z-axis
    // Apply the transform to the bounding box points
    for(i = 0; i < 3; i++) {
        for(k = 0; k < 3; k++) {
            bb_low3[i] += M[i][k]*bb_low2[k];
            bb_high3[i] += M[i][k]*bb_high2[k];
        }
    }
    // For cos theta bound, use maximum z and minimum absolute x and y
    float max_z = bb_high3[2];
    if(bb_low3[2] > max_z)
        max_z = bb_low3[2];
    // Negative maximum z means that none of the light sources are in the reflection hemisphere. Bound is 0.
    if(max_z < 0)
        return 0 * intensities;
    float min_x, min_y;
    if((bb_low3[0] < 0 && bb_high3[0] > 0) || (bb_high3[0] < 0 && bb_low3[0] > 0))
        min_x = 0;
    else
        min_x = min(fabsf(bb_low3[0]), fabsf(bb_high3[0]));
    if((bb_low3[1] < 0 && bb_high3[1] > 0) || (bb_high3[1] < 0 && bb_low3[1] > 0))
        min_y = 0;
    else
        min_y = min(fabsf(bb_low3[1]), fabsf(bb_high3[1]));
    err *= max_z / sqrtf((max_z * max_z) + (min_y * min_y) + (min_x * min_x)); // Multiply by the cosine bound
    return err * brdf * intensities; // Multiply by the constant lambertian brdf and the sum of light intensities
}

// Integrator Utility Functions
Spectrum UniformSampleAllLights(const Scene *scene,
        const Renderer *renderer, MemoryArena &arena, const Point &p,
        const Normal &n, const Vector &wo, float rayEpsilon,
        float time, BSDF *bsdf, const Sample *sample, RNG &rng,
        const LightSampleOffsets *lightSampleOffsets,
        const BSDFSampleOffsets *bsdfSampleOffsets) {
    Spectrum L(0.);
    Spectrum zero(0.);
    Vector wi; // Normal vector to sample PBRT's BRDF
    wi.x = n.x;
    wi.y = n.y;
    wi.z = n.z;
    vector<cutvertex> cut; // Heap (based on error) of vertices in the cut
    cutvertex v;
    Spectrum total; // Store the total illumination for the current cut
    v.l = scene->lighttree; // Initial cut consists only the root
    SetError(v); // Compute it's error
    Illuminate(v); // Compute it's illumination
    total = v.illumination;
    cut.push_back(v);
    while(cut.size() < 1000) {
        pop_heap(cut.begin(), cut.end()); // Get the vertex in the cut with the maximum error
        v = cut[cut.size() - 1];
        float err[3];
        float intensity[3];
        v.error.ToRGB(err);
        total.ToRGB(intensity);
        // Check if error is below the threshold. If so, we are done
        if(err[0] >= 0 && err[1] >= 0 && err[2] >= 0 && err[0] <= THRESHOLD * intensity[0] && err[1] <= THRESHOLD * intensity[1] && err[2] <= THRESHOLD * intensity[2])
            break;
        // If not, refine the cut by subdividing this vertex
        cut.pop_back(); // Remove the vertex from the cut
        total = total - v.illumination; // Remove the illumination due to this vertex
        cutvertex v1, v2; // We will add the vertex's children to the cut
        v1.l = v.l->leftChild;
        v2.l = v.l->rightChild;
        // Compute their errors and illumination
        SetError(v1);
        SetError(v2);
        Illuminate(v1);
        Illuminate(v2);
        // Add them to the cut and update the total illumination
        cut.push_back(v1);
        push_heap(cut.begin(), cut.end());
        total += v1.illumination;
        cut.push_back(v2);
        push_heap(cut.begin(), cut.end());
        total += v2.illumination;
    }
    return total;
}


Spectrum UniformSampleOneLight(const Scene *scene,
        const Renderer *renderer, MemoryArena &arena, const Point &p,
        const Normal &n, const Vector &wo, float rayEpsilon, float time,
        BSDF *bsdf, const Sample *sample, RNG &rng, int lightNumOffset,
        const LightSampleOffsets *lightSampleOffset,
        const BSDFSampleOffsets *bsdfSampleOffset) {
    // Randomly choose a single light to sample, _light_
    int nLights = int(scene->lights.size());
    if (nLights == 0) return Spectrum(0.);
    int lightNum;
    if (lightNumOffset != -1)
        lightNum = Floor2Int(sample->oneD[lightNumOffset][0] * nLights);
    else
        lightNum = Floor2Int(rng.RandomFloat() * nLights);
    lightNum = min(lightNum, nLights-1);
    Light *light = scene->lights[lightNum];

    // Initialize light and bsdf samples for single light sample
    LightSample lightSample;
    BSDFSample bsdfSample;
    if (lightSampleOffset != NULL && bsdfSampleOffset != NULL) {
        lightSample = LightSample(sample, *lightSampleOffset, 0);
        bsdfSample = BSDFSample(sample, *bsdfSampleOffset, 0);
    }
    else {
        lightSample = LightSample(rng);
        bsdfSample = BSDFSample(rng);
    }
    return (float)nLights *
        EstimateDirect(scene, renderer, arena, light, p, n, wo,
                       rayEpsilon, time, bsdf, rng, lightSample,
                       bsdfSample, BxDFType(BSDF_ALL & ~BSDF_SPECULAR));
}


Spectrum EstimateDirect(const Scene *scene, const Renderer *renderer,
        MemoryArena &arena, const Light *light, const Point &p,
        const Normal &n, const Vector &wo, float rayEpsilon, float time,
        const BSDF *bsdf, RNG &rng, const LightSample &lightSample,
        const BSDFSample &bsdfSample, BxDFType flags) {
    Spectrum Ld(0.);
    // Sample light source with multiple importance sampling
    Vector wi;
    float lightPdf, bsdfPdf;
    VisibilityTester visibility;
    Spectrum Li = light->Sample_L(p, rayEpsilon, lightSample, time,
                                  &wi, &lightPdf, &visibility);
    if (lightPdf > 0. && !Li.IsBlack()) {
        Spectrum f = bsdf->f(wo, wi, flags);
        if (!f.IsBlack() && visibility.Unoccluded(scene)) {
            // Add light's contribution to reflected radiance
            Li *= visibility.Transmittance(scene, renderer, NULL, rng, arena);
            if (light->IsDeltaLight())
                Ld += f * Li * (AbsDot(wi, n) / lightPdf);
            else {
                bsdfPdf = bsdf->Pdf(wo, wi, flags);
                float weight = PowerHeuristic(1, lightPdf, 1, bsdfPdf);
                Ld += f * Li * (AbsDot(wi, n) * weight / lightPdf);
            }
        }
    }

    // Sample BSDF with multiple importance sampling
    if (!light->IsDeltaLight()) {
        BxDFType sampledType;
        Spectrum f = bsdf->Sample_f(wo, &wi, bsdfSample, &bsdfPdf, flags,
                                    &sampledType);
        if (!f.IsBlack() && bsdfPdf > 0.) {
            float weight = 1.f;
            if (!(sampledType & BSDF_SPECULAR)) {
                lightPdf = light->Pdf(p, wi);
                if (lightPdf == 0.)
                    return Ld;
                weight = PowerHeuristic(1, bsdfPdf, 1, lightPdf);
            }
            // Add light contribution from BSDF sampling
            Intersection lightIsect;
            Spectrum Li(0.f);
            RayDifferential ray(p, wi, rayEpsilon, INFINITY, time);
            if (scene->Intersect(ray, &lightIsect)) {
                if (lightIsect.primitive->GetAreaLight() == light)
                    Li = lightIsect.Le(-wi);
            }
            else
                Li = light->Le(ray);
            if (!Li.IsBlack()) {
                Li *= renderer->Transmittance(scene, ray, NULL, rng, arena);
                Ld += f * Li * AbsDot(wi, n) * weight / bsdfPdf;
            }
        }
    }
    return Ld;
}


Spectrum SpecularReflect(const RayDifferential &ray, BSDF *bsdf,
        RNG &rng, const Intersection &isect, const Renderer *renderer,
        const Scene *scene, const Sample *sample, MemoryArena &arena) {
    Vector wo = -ray.d, wi;
    float pdf;
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    Spectrum f = bsdf->Sample_f(wo, &wi, BSDFSample(rng), &pdf,
                                BxDFType(BSDF_REFLECTION | BSDF_SPECULAR));
    Spectrum L = 0.f;
    if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, n) != 0.f) {
        // Compute ray differential _rd_ for specular reflection
        RayDifferential rd(p, wi, ray, isect.rayEpsilon);
        if (ray.hasDifferentials) {
            rd.hasDifferentials = true;
            rd.rxOrigin = p + isect.dg.dpdx;
            rd.ryOrigin = p + isect.dg.dpdy;
            // Compute differential reflected directions
            Normal dndx = bsdf->dgShading.dndu * bsdf->dgShading.dudx +
                          bsdf->dgShading.dndv * bsdf->dgShading.dvdx;
            Normal dndy = bsdf->dgShading.dndu * bsdf->dgShading.dudy +
                          bsdf->dgShading.dndv * bsdf->dgShading.dvdy;
            Vector dwodx = -ray.rxDirection - wo, dwody = -ray.ryDirection - wo;
            float dDNdx = Dot(dwodx, n) + Dot(wo, dndx);
            float dDNdy = Dot(dwody, n) + Dot(wo, dndy);
            rd.rxDirection = wi - dwodx + 2 * Vector(Dot(wo, n) * dndx +
                                                     dDNdx * n);
            rd.ryDirection = wi - dwody + 2 * Vector(Dot(wo, n) * dndy +
                                                     dDNdy * n);
        }
        PBRT_STARTED_SPECULAR_REFLECTION_RAY(const_cast<RayDifferential *>(&rd));
        Spectrum Li = renderer->Li(scene, rd, sample, rng, arena);
        L = f * Li * AbsDot(wi, n) / pdf;
        PBRT_FINISHED_SPECULAR_REFLECTION_RAY(const_cast<RayDifferential *>(&rd));
    }
    return L;
}


Spectrum SpecularTransmit(const RayDifferential &ray, BSDF *bsdf,
        RNG &rng, const Intersection &isect, const Renderer *renderer,
        const Scene *scene, const Sample *sample, MemoryArena &arena) {
    Vector wo = -ray.d, wi;
    float pdf;
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    Spectrum f = bsdf->Sample_f(wo, &wi, BSDFSample(rng), &pdf,
                               BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR));
    Spectrum L = 0.f;
    if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, n) != 0.f) {
        // Compute ray differential _rd_ for specular transmission
        RayDifferential rd(p, wi, ray, isect.rayEpsilon);
        if (ray.hasDifferentials) {
            rd.hasDifferentials = true;
            rd.rxOrigin = p + isect.dg.dpdx;
            rd.ryOrigin = p + isect.dg.dpdy;
        
            float eta = bsdf->eta;
            Vector w = -wo;
            if (Dot(wo, n) < 0) eta = 1.f / eta;
        
            Normal dndx = bsdf->dgShading.dndu * bsdf->dgShading.dudx + bsdf->dgShading.dndv * bsdf->dgShading.dvdx;
            Normal dndy = bsdf->dgShading.dndu * bsdf->dgShading.dudy + bsdf->dgShading.dndv * bsdf->dgShading.dvdy;
        
            Vector dwodx = -ray.rxDirection - wo, dwody = -ray.ryDirection - wo;
            float dDNdx = Dot(dwodx, n) + Dot(wo, dndx);
            float dDNdy = Dot(dwody, n) + Dot(wo, dndy);
        
            float mu = eta * Dot(w, n) - Dot(wi, n);
            float dmudx = (eta - (eta*eta*Dot(w,n))/Dot(wi, n)) * dDNdx;
            float dmudy = (eta - (eta*eta*Dot(w,n))/Dot(wi, n)) * dDNdy;
        
            rd.rxDirection = wi + eta * dwodx - Vector(mu * dndx + dmudx * n);
            rd.ryDirection = wi + eta * dwody - Vector(mu * dndy + dmudy * n);
        }
        PBRT_STARTED_SPECULAR_REFRACTION_RAY(const_cast<RayDifferential *>(&rd));
        Spectrum Li = renderer->Li(scene, rd, sample, rng, arena);
        L = f * Li * AbsDot(wi, n) / pdf;
        PBRT_FINISHED_SPECULAR_REFRACTION_RAY(const_cast<RayDifferential *>(&rd));
    }
    return L;
}


Distribution1D *ComputeLightSamplingCDF(const Scene *scene) {
    uint32_t nLights = int(scene->lights.size());
    Assert(nLights > 0);
    vector<float>lightPower(nLights, 0.f);
    for (uint32_t i = 0; i < nLights; ++i)
        lightPower[i] = scene->lights[i]->Power(scene).y();
    return new Distribution1D(&lightPower[0], nLights);
}


