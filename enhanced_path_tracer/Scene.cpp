//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                pdf = 1.0f / emit_area_sum;  //Set the PDF to the reciprocal of the joint light source area.
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    if (depth > this->maxDepth) {
        return Vector3f();
    }

    Vector3f wo = -ray.direction;
    Intersection pos = Scene::intersect(ray);
    Vector3f p = pos.coords;
    Vector3f N = pos.normal;

    //Miss hitting objects test 
    if (!pos.happened)
        return Vector3f();

    //Calculate the light source illumination

    Vector3f L_lig = Vector3f();
    if (pos.m->hasEmission()) {
        L_lig = pos.m->getEmission() * dotProduct(wo, pos.normal);
    }

    //Calculate the direct illumination
    Vector3f L_dir = Vector3f();
    Intersection inter = Intersection();
    float pdf_light = 0;

    Scene::sampleLight(inter, pdf_light);

    Material* ma = pos.m;
    Vector3f x = inter.coords;
    Vector3f NN = inter.normal;
    Vector3f emit = inter.emit;
    Vector3f ws = (x - p).normalized();
    float dis2 = std::pow((x - p).norm(), 2);

    //Light source visibility test
    Vector3f testOri = p + N * EPSILON;
    Ray testRay = Ray(testOri, ws);
    Intersection testInter = Scene::intersect(testRay);
    float testDis = (testInter.coords - x).norm();

    if (testInter.happened && testInter.m->hasEmission())
    {
        L_dir = emit * ma->eval(ws, wo, N) * dotProduct(ws, N) * dotProduct(-ws, NN) / (dis2 * pdf_light);
    }

    //Calculate indirect illumination
    Vector3f L_indir = Vector3f();

    //Test Russian Roulette
    float pro = get_random_float();
    if (pro < RussianRoulette) {
        Vector3f wi = (ma->sample(wo, N)).normalized();
        Intersection indirInter = Scene::intersect(Ray(p, wi));

        //If the ray hits an object that does not emit light
        if (indirInter.happened && !indirInter.m->hasEmission())
        {
            //debug message
            //float pdf_val = ma->pdf(wi, wo, N);
            //Vector3f eval_val = ma->eval(wi, wo, N);
            //float cos_term = dotProduct(wi, N);
            //Vector3f weight = eval_val * cos_term / pdf_val;
            //if (weight.x > 10.0f || weight.y > 10.0f || weight.z > 10.0f) {
            //    std::cout << "Huge weight: " << weight << " pdf=" << pdf_val << " eval=" << eval_val << std::endl;
            //}

            L_indir = castRay(Ray(p, wi), depth + 1) * ma->eval(wi, wo, N) * dotProduct(wi, N) / (ma->pdf(wi, wo, N) * RussianRoulette);
        }
		//std::cout << "L_indir: " << L_indir << "\n";
        //std::cout<<"L_dir eval="<<ma->eval(wi, wo, N)<<std::endl;
    }

    return L_lig + L_dir + L_indir;
    //Custom implement ends
}