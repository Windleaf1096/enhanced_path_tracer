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
    float testDis = (testInter.coords - x).norm();  //The testDis variable is not used.


	//L_dir dubug message

    static int debug_count = 0;
    const int max_debug_print = 20; 

    //Check if the hit point is emissive to determine the light source.
    if (testInter.happened && testInter.m->hasEmission())
    {
        //Vector3f eval_val = ma->eval(ws, wo, N);
        //float cos_ws_N = dotProduct(ws, N);
        //float cos_negws_NN = dotProduct(-ws, NN);
        //float pdf_light_val = pdf_light;

        //Vector3f L_dir_raw = emit * eval_val * cos_ws_N * cos_negws_NN / (dis2 * pdf_light_val);

        //if (debug_count < max_debug_print)
        //{
        //    printf("\n[DirectLight] Sample #%d\n", debug_count++);
        //    printf("  emit = (%.4f, %.4f, %.4f)\n", emit.x, emit.y, emit.z);
        //    printf("  eval = (%.6f, %.6f, %.6f)\n", eval_val.x, eval_val.y, eval_val.z);
        //    printf("  ws = (%.4f, %.4f, %.4f), N = (%.4f, %.4f, %.4f)\n", ws.x, ws.y, ws.z, N.x, N.y, N.z);
        //    printf("  dot(ws, N) = %.6f\n", cos_ws_N);
        //    printf("  -ws = (%.4f, %.4f, %.4f), NN = (%.4f, %.4f, %.4f)\n", -ws.x, -ws.y, -ws.z, NN.x, NN.y, NN.z);
        //    printf("  dot(-ws, NN) = %.6f\n", cos_negws_NN);
        //    printf("  dis2 = %.6f, pdf_light = %.6f\n", dis2, pdf_light_val);
        //    printf("  L_dir = (%.6f, %.6f, %.6f)\n", L_dir_raw.x, L_dir_raw.y, L_dir_raw.z);

        //    if (L_dir_raw.x > 10.0f || L_dir_raw.y > 10.0f || L_dir_raw.z > 10.0f)
        //        printf("  [WARNING] L_dir component exceeds 10!\n");
        //}

        L_dir = emit * ma->eval(ws, wo, N) * dotProduct(ws, N) * dotProduct(-ws, NN) / (dis2 * pdf_light);
    }

    // Calculate indirect illumination
    Vector3f L_indir = Vector3f();

    // Test Russian Roulette
    float pro = get_random_float();
    if (pro < RussianRoulette) {
        Vector3f wi = (ma->sample(wo, N)).normalized();
        Intersection indirInter = Scene::intersect(Ray(p, wi));

        //If the ray hits an object that does not emit light
        if (indirInter.happened && !indirInter.m->hasEmission())
        {
            L_indir = castRay(Ray(p, wi), depth + 1) * ma->eval(wi, wo, N) * dotProduct(wi, N) / (ma->pdf(wi, wo, N) * RussianRoulette);
        }
    }

    return L_lig + L_dir + L_indir;
    //Custom implement ends
}