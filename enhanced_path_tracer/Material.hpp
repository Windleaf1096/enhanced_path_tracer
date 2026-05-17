//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

#include "Vector.hpp"

enum MaterialType { DIFFUSE,MICROFACET,GGX};

class Material{
private:

    // Compute reflection direction
    Vector3f reflect(const Vector3f &I, const Vector3f &N) const
    {
        return I - 2 * dotProduct(I, N) * N;
    }

    // Compute refraction direction using Snell's law
    //
    // We need to handle with care the two possible situations:
    //
    //    - When the ray is inside the object
    //
    //    - When the ray is outside.
    //
    // If the ray is outside, you need to make cosi positive cosi = -N.I
    //
    // If the ray is inside, you need to invert the refractive indices and negate the normal N
    Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        Vector3f n = N;
        if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
    }

    // Compute Fresnel equation
    //
    // \param I is the incident view direction
    //
    // \param N is the normal at the intersection point
    //
    // \param ior is the material refractive index
    //
    // \param[out] kr is the amount of light reflected
    void fresnel(const Vector3f &I, const Vector3f &N, const float &ior, float &kr) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        if (cosi > 0) {  std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1) {
            kr = 1;
        }
        else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            kr = (Rs * Rs + Rp * Rp) / 2;
        }
        // As a consequence of the conservation of energy, transmittance is given by:
        // kt = 1 - kr;
    }

    Vector3f toWorld(const Vector3f &a, const Vector3f &N){
        Vector3f B, C;
        if (std::fabs(N.x) > std::fabs(N.y)){
            float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
            C = Vector3f(N.z * invLen, 0.0f, -N.x *invLen);
        }
        else {
            float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
            C = Vector3f(0.0f, N.z * invLen, -N.y *invLen);
        }
        B = crossProduct(C, N);
        return a.x * B + a.y * C + a.z * N;
    }

    // Sample a halfvector according to GGX distribution
    Vector3f sampleGGX(const Vector3f& N, float alpha) {
        float r1 = get_random_float();
        float r2 = get_random_float();

        // theta of halfvector
        float theta = std::atan(alpha * std::sqrt(r1) / std::sqrt(1.0f - r1));
        float phi = 2.0f * M_PI * r2;

        float sinTheta = std::sin(theta);
        float cosTheta = std::cos(theta);

        Vector3f localH(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
        return toWorld(localH, N);
    }

	// GGX NDF function
	float D_GGX(const Vector3f& N, const Vector3f& H, float alpha) {
		float cosThetaH = dotProduct(N, H);
		if (cosThetaH <= 0.0f) return 0.0f;
		float alpha2 = alpha * alpha;
		float denom = cosThetaH * cosThetaH * (alpha2 - 1.0f) + 1.0f;
		return alpha2 / (M_PI * denom * denom);
	}



public:
    MaterialType m_type;
    //Vector3f m_color;
    Vector3f m_emission;
    float ior;
    Vector3f Kd, Ks;
    float specularExponent;
    //Texture tex;
	float roughness; // roughness parameter for GGX
	float alpha; // roughness parameter for GGX

    inline Material(MaterialType t = DIFFUSE, Vector3f e = Vector3f(0, 0, 0), float r = 0.2f);
    inline MaterialType getType();
    //inline Vector3f getColor();
    inline Vector3f getColorAt(double u, double v);
    inline Vector3f getEmission();
    inline bool hasEmission();

    // sample a ray by Material properties
    inline Vector3f sample(const Vector3f &wi, const Vector3f &N);
    // given a ray, calculate the PdF of this ray
    inline float pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);
    // given a ray, calculate the contribution of this ray
    inline Vector3f eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);

};

Material::Material(MaterialType t, Vector3f e,float r){
    m_type = t;
    m_emission = e;
    Kd = Vector3f(0.0f);
    Ks = Vector3f(0.0f);
    ior = 1.3f;
    specularExponent = 50.0f;
    roughness = r;
    alpha = roughness * roughness;
}

MaterialType Material::getType(){return m_type;}
///Vector3f Material::getColor(){return m_color;}
Vector3f Material::getEmission() {return m_emission;}
bool Material::hasEmission() {

    if (m_emission.norm() > EPSILON) {
        return true;
    }
    else 
        return false;
}

Vector3f Material::getColorAt(double u, double v) {
    return Vector3f();
}

//Given an incoming ray and the normal at the intersection point, sample a ray according to the material properties
Vector3f Material::sample(const Vector3f &wi, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1);
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r*std::cos(phi), r*std::sin(phi), z);
            return toWorld(localRay, N);
            
            break;
        }

        case MICROFACET:
        {
            // uniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1);
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r * std::cos(phi), r * std::sin(phi), z);
            return toWorld(localRay, N);

            break;
        }

        //GGX NDF-based importance sampling
        case GGX: {
            // Sample half-vector h according to GGX NDF
			Vector3f h = sampleGGX(N, alpha);  

            //Takes the incident direction toward the surface
			Vector3f wo = reflect(wi, h);  

            //Returns the outgoing direction leaving the surface
            return -wo;
        }

    }
}

//given a ray, calculate the pdf of this ray
float Material::pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample probability 1 / (2 * PI)
            if (dotProduct(wi, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }

        //MICROFACET uses the same pdf function
        case MICROFACET:
        {
            // uniform sample probability 1 / (2 * PI)
            if (dotProduct(wi, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }

		//GGX 
        case GGX: {
			//Calculate the half-vector h
			Vector3f h = normalize(wi + wo);
			float D = D_GGX(N, h, alpha);
            float cosThetaH = std::max(0.0f, dotProduct(N, h));
			float cosOH = std::max(0.0f, dotProduct(wo, h));

			if (cosOH < 1e-6f) return 0.0f; // Avoid division by zero

			//PDF for sampling wo given wi
			float pdf_wi = D * cosThetaH / (4.0f * cosOH);
			return pdf_wi;
        }

    }
}


//calculate the BRDF
Vector3f Material::eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // calculate the contribution of diffuse   model
            float cosalpha = dotProduct(N, wo);
            //float cosalpha = dotProduct(N, wi);
            if (cosalpha > 0.0f) {
                Vector3f diffuse = Kd / M_PI;
                return diffuse;
            }
            else
                return Vector3f(0.0f);
            break;
        }

        //eval BRDF for MICROFACET materials
        case MICROFACET:
        {
            float cos_wo = dotProduct(N, wo);
            float cos_wi = dotProduct(N, wi);
            if (cos_wi <= 0.0f || cos_wo <= 0.0f)
                return Vector3f(0.0f);

            Vector3f H_unnormalized = wi + wo;
            if (H_unnormalized.norm() < 1e-6f)
                return Kd / M_PI;  
            Vector3f H = normalize(wi + wo);

            float alpha = std::sqrt(2.0f / (specularExponent + 2.0f));

            float cos_theta_h = dotProduct(N, H);
            if (cos_theta_h <= 0.0f) return Vector3f(0.0f);
            float alpha2 = alpha * alpha;
            float denomD = cos_theta_h * cos_theta_h * (alpha2 - 1.0f) + 1.0f;
            float D = alpha2 / (M_PI * denomD * denomD);

            auto G1 = [&](const Vector3f& v, const Vector3f& n) -> float {
                float cos_v = dotProduct(v, n);
                if (cos_v <= 0.0f) return 0.0f;
                float tan2 = (1.0f - cos_v * cos_v) / (cos_v * cos_v);
                return 2.0f / (1.0f + std::sqrt(1.0f + alpha2 * tan2));
                };
            float G = G1(wi, N) * G1(wo, N);

            float F = 0.0f;
            fresnel(wi, H, ior, F);

            float denom = 4.0f * cos_wi * cos_wo;
            if (denom <= 0.0f) return Vector3f(0.0f);
            Vector3f specular = Ks * (D * G * F / denom);

            Vector3f diffuse = Kd / M_PI;

            return diffuse + specular;
            break;
        }

		//eval BRDF for GGX materials
        case GGX: {
            // GGX BRDF evaluation
        }
    }
}

#endif //RAYTRACING_MATERIAL_H
