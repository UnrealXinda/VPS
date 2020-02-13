// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "ProceduralMeshComponent.h"
#include "ProceduralClothComponent.generated.h"

struct FClothParticle
{
	int     bFree;
	FVector OldLoc;
	FVector NewLoc;
};

/**
 * 
 */
UCLASS(hidecategories = (Object, LOD), editinlinenew, meta = (BlueprintSpawnableComponent), ClassGroup = Rendering, DisplayName = "ProceduralClothComponent")
class VPS_API UProceduralClothComponent : public UProceduralMeshComponent
{
	GENERATED_BODY()

public:

	UPROPERTY(EditAnywhere, BlueprintReadWrite, meta = (ClampMin = 2), Category = "Cloth Geometry")
	int32 HorizontalVertexCount;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, meta = (ClampMin = 2), Category = "Cloth Geometry")
	int32 VerticalVertexCount;
	
	UPROPERTY(EditAnywhere, BlueprintReadWrite, meta = (ClampMin = 0), Category = "Cloth Geometry")
	float HorizontalDistance;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, meta = (ClampMin = 0), Category = "Cloth Geometry")
	float VerticalDistance;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, meta = (ClampMin = 0.0, ClampMax = 1.0), Category = "Cloth Physics")
	float ClothGravityScale;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, meta = (ClampMin = 0.005, ClampMax = 0.1), Category = "Cloth Physics")
	float SubstepTime;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, meta = (ClampMin = 1), Category = "Cloth Physics")
	int32 SolverIterationCount;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Cloth Physics")
	uint8 bUseGPUAcceleration : 1;

public:

	UProceduralClothComponent(const FObjectInitializer& ObjectInitializer);

	UFUNCTION(BlueprintCallable)
	void InitClothComponent();

	void InitClothParticles();
	void UpdateProceduralMesh(bool bVertexOnly);
	void TickClothComponent(float DeltaTime);

	// CPU parallel solver
	void PerformSubstep(float InSubstepTime, const FVector& Gravity);
	void VerletIntegrate(float InSubstepTime, const FVector& Gravity);
	void SolveConstraints();

	// GPU parallel solver
	void PerformSubstepParallel(float InSubstepTime, const FVector& Gravity);
	void VerletIntegrateParallel(float InSubstepTime, const FVector& Gravity);
	void SolveConstraintsParallel();
	void ReleaseBufferResources();

	virtual void OnRegister() override;
	virtual void TickComponent(float DeltaTime, enum ELevelTick TickType, FActorComponentTickFunction *ThisTickFunction) override;

protected:

	/** Amount of time 'left over' from last tick */
	float TimeRemainder;

	TResourceArray<FClothParticle> ClothParticles;  // Declared as TResourceArray to be compatible with GPU solver

	FStructuredBufferRHIRef    ParticlesStructuredBuffer;
	FUnorderedAccessViewRHIRef ParticlesStructuredBufferUAV;

};
