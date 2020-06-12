// Fill out your copyright notice in the Description page of Project Settings.


#include "ProceduralClothComponent.h"
#include "Runtime/Core/Public/Async/ParallelFor.h"
#include "RenderCore/Public/GlobalShader.h"
#include "RenderCore/Public/ShaderParameterUtils.h"
#include "RenderCore/Public/ShaderParameterMacros.h"

#include "Public/GlobalShader.h"
#include "Public/PipelineStateCache.h"
#include "Public/RHIStaticStates.h"
#include "Public/SceneUtils.h"
#include "Public/SceneInterface.h"
#include "Public/ShaderParameterUtils.h"
#include "Public/Logging/MessageLog.h"
#include "Public/Internationalization/Internationalization.h"
#include "Public/StaticBoundShaderState.h"
#include "RHI/Public/RHICommandList.h"

BEGIN_GLOBAL_SHADER_PARAMETER_STRUCT(FClothVerletComputeShaderParameters, )
	SHADER_PARAMETER(int,     SizeX)
	SHADER_PARAMETER(int,     SizeY)
	SHADER_PARAMETER(float,   SubStepSqr)
	SHADER_PARAMETER(FVector, Gravity)
END_GLOBAL_SHADER_PARAMETER_STRUCT()
IMPLEMENT_GLOBAL_SHADER_PARAMETER_STRUCT(FClothVerletComputeShaderParameters, "ClothVerletUniform");

BEGIN_GLOBAL_SHADER_PARAMETER_STRUCT(FClothConstraintComputeShaderParameters, )
	SHADER_PARAMETER(int,     SizeX)
	SHADER_PARAMETER(int,     SizeY)
	SHADER_PARAMETER(int,     Direction)   // 0 for solving horizontal constraint. 1 for solving vertical constraint
	SHADER_PARAMETER(float,   DesiredDistance)
	SHADER_PARAMETER(float,   GroundZ)
	SHADER_PARAMETER(float,   SphereRadius)
	SHADER_PARAMETER(FVector, SphereLoc)
END_GLOBAL_SHADER_PARAMETER_STRUCT()
IMPLEMENT_GLOBAL_SHADER_PARAMETER_STRUCT(FClothConstraintComputeShaderParameters, "ClothConstraintUniform");

BEGIN_GLOBAL_SHADER_PARAMETER_STRUCT(FClothNormalComputeShaderParameters, )
	SHADER_PARAMETER(int, SizeX)
	SHADER_PARAMETER(int, SizeY)
END_GLOBAL_SHADER_PARAMETER_STRUCT()
IMPLEMENT_GLOBAL_SHADER_PARAMETER_STRUCT(FClothNormalComputeShaderParameters, "ClothNormalUniform");

class FClothVerletComputeShader : public FGlobalShader
{

	DECLARE_SHADER_TYPE(FClothVerletComputeShader, Global)

public:
	FClothVerletComputeShader() {}
	FClothVerletComputeShader(const ShaderMetaType::CompiledShaderInitializerType& Initializer)
		: FGlobalShader(Initializer)
	{
		OutputParticleBuffer.Bind(Initializer.ParameterMap, TEXT("OutputParticleBuffer"));
	}

	static bool ShouldCompilePermutation(const FGlobalShaderPermutationParameters& Parameters)
	{
		return IsFeatureLevelSupported(Parameters.Platform, ERHIFeatureLevel::SM5);
	}

	static bool ShouldCache(EShaderPlatform Platform)
	{
		return IsFeatureLevelSupported(Platform, ERHIFeatureLevel::SM5);
	}

	void BindShaderBuffers(FRHICommandList& RHICmdList, FUnorderedAccessViewRHIRef OutputParticleBufferUAV)
	{
		FRHIComputeShader* ComputeShaderRHI = RHICmdList.GetBoundComputeShader();
		SetUAVParameter(RHICmdList, ComputeShaderRHI, OutputParticleBuffer, OutputParticleBufferUAV);
	}

	void UnbindShaderBuffers(FRHICommandList& RHICmdList)
	{
		FRHIComputeShader* ComputeShaderRHI = RHICmdList.GetBoundComputeShader();
		SetUAVParameter(RHICmdList, ComputeShaderRHI, OutputParticleBuffer, FUnorderedAccessViewRHIRef());
	}

	void SetShaderParameters(FRHICommandList& RHICmdList, const FClothVerletComputeShaderParameters& Parameters)
	{
		FRHIComputeShader* ComputeShaderRHI = RHICmdList.GetBoundComputeShader();
		SetUniformBufferParameterImmediate(RHICmdList, ComputeShaderRHI, GetUniformBufferParameter<FClothVerletComputeShaderParameters>(), Parameters);
	}

private:

	LAYOUT_FIELD(FShaderResourceParameter, OutputParticleBuffer);
};

class FClothConstraintComputeShader : public FGlobalShader
{

	DECLARE_SHADER_TYPE(FClothConstraintComputeShader, Global)

public:
	FClothConstraintComputeShader() {}
	FClothConstraintComputeShader(const ShaderMetaType::CompiledShaderInitializerType& Initializer)
		: FGlobalShader(Initializer)
	{
		OutputParticleBuffer.Bind(Initializer.ParameterMap, TEXT("OutputParticleBuffer"));
		OutputPositionBuffer.Bind(Initializer.ParameterMap, TEXT("OutputPositionBuffer"));
	}

	static bool ShouldCompilePermutation(const FGlobalShaderPermutationParameters& Parameters)
	{
		return IsFeatureLevelSupported(Parameters.Platform, ERHIFeatureLevel::SM5);
	}

	static bool ShouldCache(EShaderPlatform Platform)
	{
		return IsFeatureLevelSupported(Platform, ERHIFeatureLevel::SM5);
	}

	void BindShaderBuffers(FRHICommandList& RHICmdList, FUnorderedAccessViewRHIRef OutputParticleBufferUAV, FUnorderedAccessViewRHIRef OutputPositionBufferUAV)
	{
		FRHIComputeShader* ComputeShaderRHI = RHICmdList.GetBoundComputeShader();
		SetUAVParameter(RHICmdList, ComputeShaderRHI, OutputParticleBuffer, OutputParticleBufferUAV);
		SetUAVParameter(RHICmdList, ComputeShaderRHI, OutputPositionBuffer, OutputPositionBufferUAV);
	}

	void UnbindShaderBuffers(FRHICommandList& RHICmdList)
	{
		FRHIComputeShader* ComputeShaderRHI = RHICmdList.GetBoundComputeShader();
		SetUAVParameter(RHICmdList, ComputeShaderRHI, OutputParticleBuffer, FUnorderedAccessViewRHIRef());
		SetUAVParameter(RHICmdList, ComputeShaderRHI, OutputPositionBuffer, FUnorderedAccessViewRHIRef());
	}

	void SetShaderParameters(FRHICommandList& RHICmdList, const FClothConstraintComputeShaderParameters& Parameters)
	{
		FRHIComputeShader* ComputeShaderRHI = RHICmdList.GetBoundComputeShader();
		SetUniformBufferParameterImmediate(RHICmdList, ComputeShaderRHI, GetUniformBufferParameter<FClothConstraintComputeShaderParameters>(), Parameters);
	}

private:

	LAYOUT_FIELD(FShaderResourceParameter, OutputParticleBuffer);
	LAYOUT_FIELD(FShaderResourceParameter, OutputPositionBuffer);
};

class FClothNormalComputeShader : public FGlobalShader
{

	DECLARE_SHADER_TYPE(FClothNormalComputeShader, Global)

public:
	FClothNormalComputeShader() {}
	FClothNormalComputeShader(const ShaderMetaType::CompiledShaderInitializerType& Initializer)
		: FGlobalShader(Initializer)
	{
		OutputPositionBuffer.Bind(Initializer.ParameterMap, TEXT("OutputPositionBuffer"));
		OutputNormalBuffer.Bind(Initializer.ParameterMap, TEXT("OutputNormalBuffer"));
	}

	static bool ShouldCompilePermutation(const FGlobalShaderPermutationParameters& Parameters)
	{
		return IsFeatureLevelSupported(Parameters.Platform, ERHIFeatureLevel::SM5);
	}

	static bool ShouldCache(EShaderPlatform Platform)
	{
		return IsFeatureLevelSupported(Platform, ERHIFeatureLevel::SM5);
	}

	void BindShaderBuffers(FRHICommandList& RHICmdList, FUnorderedAccessViewRHIRef OutputPositionBufferUAV, FUnorderedAccessViewRHIRef OutputNormalBufferUAV)
	{
		FRHIComputeShader* ComputeShaderRHI = RHICmdList.GetBoundComputeShader();
		SetUAVParameter(RHICmdList, ComputeShaderRHI, OutputPositionBuffer, OutputPositionBufferUAV);
		SetUAVParameter(RHICmdList, ComputeShaderRHI, OutputNormalBuffer,   OutputNormalBufferUAV);
	}

	void UnbindShaderBuffers(FRHICommandList& RHICmdList)
	{
		FRHIComputeShader* ComputeShaderRHI = RHICmdList.GetBoundComputeShader();
		SetUAVParameter(RHICmdList, ComputeShaderRHI, OutputPositionBuffer, FUnorderedAccessViewRHIRef());
		SetUAVParameter(RHICmdList, ComputeShaderRHI, OutputNormalBuffer,   FUnorderedAccessViewRHIRef());
	}

	void SetShaderParameters(FRHICommandList& RHICmdList, const FClothNormalComputeShaderParameters& Parameters)
	{
		FRHIComputeShader* ComputeShaderRHI = RHICmdList.GetBoundComputeShader();
		SetUniformBufferParameterImmediate(RHICmdList, ComputeShaderRHI, GetUniformBufferParameter<FClothNormalComputeShaderParameters>(), Parameters);
	}

private:

	LAYOUT_FIELD(FShaderResourceParameter, OutputPositionBuffer);
	LAYOUT_FIELD(FShaderResourceParameter, OutputNormalBuffer);
};

IMPLEMENT_GLOBAL_SHADER(FClothVerletComputeShader,     "/Plugin/VPS/ClothConstraintComputeShader.usf", "ComputeClothVerlet",   SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FClothNormalComputeShader,     "/Plugin/VPS/ClothConstraintComputeShader.usf", "ComputeClothNormal",   SF_Compute);
IMPLEMENT_GLOBAL_SHADER(FClothConstraintComputeShader, "/Plugin/VPS/ClothConstraintComputeShader.usf", "SolveClothConstraint", SF_Compute);

namespace
{
	template <class T>
	FORCEINLINE void ReleaseBufferResource(T Resource)
	{
		if (Resource.IsValid())
		{
			Resource->Release();
		}
	}

	void SolveDistanceConstraint(FClothParticle& ParticleA, FClothParticle& ParticleB, float DesiredDistance)
	{
		// Find current vector between particles
		FVector Delta = ParticleB.NewLoc - ParticleA.NewLoc;

		float CurrentDistance = Delta.Size();
		float ErrorFactor = (CurrentDistance - DesiredDistance) / CurrentDistance;

		// Only move free particles to satisfy constraints
		if (ParticleA.bFree && ParticleB.bFree)
		{
			ParticleA.NewLoc += ErrorFactor * 0.5f * Delta;
			ParticleB.NewLoc -= ErrorFactor * 0.5f * Delta;
		}
		else if (ParticleA.bFree)
		{
			ParticleA.NewLoc += ErrorFactor * Delta;
		}
		else if (ParticleB.bFree)
		{
			ParticleB.NewLoc -= ErrorFactor * Delta;
		}
	}

#if VPS_CLOTH_DEBUG

	void SolveSphericalConstraint(FClothParticle& Particle, float SphereRadius, FVector SphereLoc)
	{
		FVector Dir = Particle.NewLoc - SphereLoc;
		float DistSqr = FVector::DotProduct(Dir, Dir);

		if (DistSqr < SphereRadius * SphereRadius)
		{
			Particle.NewLoc = SphereLoc + SphereRadius * Dir.GetSafeNormal();
		}
	}

	inline void SolveGroundConstraint(FClothParticle& Particle, float GroundZ)
	{
		Particle.NewLoc.Z = FMath::Max(GroundZ, Particle.NewLoc.Z);
	}

#endif
}

UProceduralClothComponent::UProceduralClothComponent(const FObjectInitializer& ObjectInitializer) :
	Super(ObjectInitializer)
{
	HorizontalVertexCount = 30;
	VerticalVertexCount = 30;

	HorizontalDistance = 2.f;
	VerticalDistance = 3.f;

	ClothGravityScale = 1.0f;
	SubstepTime = 0.02f;
	SolverIterationCount = 1;

	GroundZ = -10000.0f;
	SphereRadius = 50.0f;
	SphereLoc = FVector(0, 0, -10000);

	PrimaryComponentTick.bCanEverTick = true;
	PrimaryComponentTick.bStartWithTickEnabled = true;

	bTickInEditor = true;
	bAutoActivate = true;
}

void UProceduralClothComponent::OnRegister()
{
	Super::OnRegister();

	InitClothComponent();
}

void UProceduralClothComponent::InitClothComponent()
{
	InitClothParticles();
	UpdateProceduralMesh(false);

	if (bUseGPUAcceleration)
	{
		// Release buffer before allocating new ones
		ReleaseBufferResources();

		const int32 VertexCount = HorizontalVertexCount * VerticalVertexCount;

		TResourceArray<FClothParticle> TempResourceArray;
		TempResourceArray.AddUninitialized(VertexCount);
		FMemory::Memcpy(TempResourceArray.GetData(), ClothParticles.GetData(), sizeof(FClothParticle) * VertexCount);

		FRHIResourceCreateInfo CreateInfo(&TempResourceArray);
		ParticlesStructuredBuffer = RHICreateStructuredBuffer(
			sizeof(FClothParticle),                   // Stride
			sizeof(FClothParticle) * VertexCount,     // Size
			BUF_UnorderedAccess,                      // Usage
			CreateInfo                                // Create info
		);
		ParticlesStructuredBufferUAV = RHICreateUnorderedAccessView(ParticlesStructuredBuffer, true, false);

		CreateInfo.ResourceArray = nullptr;
		PositionsStructuredBuffer = RHICreateStructuredBuffer(
			sizeof(FVector),                          // Stride
			sizeof(FVector) * VertexCount,            // Size
			BUF_UnorderedAccess,                      // Usage
			CreateInfo                                // Create info
		);
		PositionsStructuredBufferUAV = RHICreateUnorderedAccessView(PositionsStructuredBuffer, true, false);

		NormalsStructuredBuffer = RHICreateStructuredBuffer(
			sizeof(FVector),                          // Stride
			sizeof(FVector) * VertexCount,            // Size
			BUF_UnorderedAccess,                      // Usage
			CreateInfo                                // Create info
		);
		NormalsStructuredBufferUAV = RHICreateUnorderedAccessView(NormalsStructuredBuffer, true, false);
	}
}

void UProceduralClothComponent::InitClothParticles()
{
	const int32 VertexCount = HorizontalVertexCount * VerticalVertexCount;

	ClothParticles.Reset();
	ClothPositions.Reset();
	ClothNormals.Reset();
	ClothParticles.AddUninitialized(VertexCount);
	ClothPositions.AddUninitialized(VertexCount);
	ClothNormals.AddUninitialized(VertexCount);

	for (int32 Y = 0; Y < VerticalVertexCount; ++Y)
	{
		for (int32 X = 0; X < HorizontalVertexCount; ++X)
		{
			float LocX = X * HorizontalDistance;
			float LocY = Y * VerticalDistance;
			float LocZ = 0.0f;
			FVector Loc(LocX, LocY, LocZ);

			int32 Index = Y * HorizontalVertexCount + X;

			FClothParticle& ClothParticle = ClothParticles[Index];

			ClothParticle.OldLoc = Loc;
			ClothParticle.NewLoc = Loc;
			ClothParticle.bFree = !FixedParticleIndices.Contains(Index);
			//ClothParticle.bFree = (X - HorizontalVertexCount + 1) != 0 && (Y - VerticalVertexCount + 1) != 0;

			ClothPositions[Index] = Loc;
			ClothNormals[Index] = FVector::UpVector;
		}
	}
}

void UProceduralClothComponent::UpdateProceduralMesh(bool bPositionAndNormalOnly)
{
	const int32 VertexCount = HorizontalVertexCount * VerticalVertexCount;

	TArray<FVector2D> UVs;
	TArray<FColor> VertexColors;
	TArray<FProcMeshTangent> Tangents;

	if (!bPositionAndNormalOnly)
	{
		const int32 TriangleCount = (HorizontalVertexCount - 1) * (VerticalVertexCount - 1) * 6;		
		TArray<int32> Triangles;

		UVs.AddUninitialized(VertexCount);
		Triangles.AddUninitialized(TriangleCount);

		for (int32 Y = 0; Y < VerticalVertexCount; ++Y)
		{
			for (int32 X = 0; X < HorizontalVertexCount; ++X)
			{
				float U = StaticCast<float>(X) / HorizontalVertexCount - 1;
				float V = StaticCast<float>(Y) / VerticalVertexCount - 1;
				int32 Index = Y * HorizontalVertexCount + X;

				UVs[Index] = FVector2D(U, V);
			}
		}

		for (int32 Y = 0, Index = 0; Y < VerticalVertexCount - 1; ++Y)
		{
			for (int32 X = 0; X < HorizontalVertexCount - 1; ++X)
			{
				int A = Y * HorizontalVertexCount + X;
				int B = A + HorizontalVertexCount;
				int C = A + HorizontalVertexCount + 1;
				int D = A + 1;

				Triangles[Index++] = A;
				Triangles[Index++] = B;
				Triangles[Index++] = C;

				Triangles[Index++] = A;
				Triangles[Index++] = C;
				Triangles[Index++] = D;
			}
		}

		ClearMeshSection(0);
		CreateMeshSection(0, ClothPositions, Triangles, ClothNormals, UVs, VertexColors, Tangents, false);
	}

	else
	{
		UpdateMeshSection(0, ClothPositions, ClothNormals, UVs, VertexColors, Tangents);
	}
}

void UProceduralClothComponent::PerformSubstep(float InSubstepTime, const FVector& Gravity)
{
	VerletIntegrate(InSubstepTime, Gravity);
	SolveConstraints();
	ComputeNormals();
}

void UProceduralClothComponent::VerletIntegrate(float InSubstepTime, const FVector& Gravity)
{
	const float SubstepTimeSqr = InSubstepTime * InSubstepTime;

	ParallelFor(ClothParticles.Num(), [&](int32 Index)
	{
		FClothParticle& Particle = ClothParticles[Index];

		if (Particle.bFree)
		{
			const FVector DeltaLoc = Particle.NewLoc - Particle.OldLoc;
			const FVector NewLoc = Particle.NewLoc + DeltaLoc + (SubstepTimeSqr * Gravity);

			Particle.OldLoc = Particle.NewLoc;
			Particle.NewLoc = NewLoc;
		}
	});
}

void UProceduralClothComponent::SolveConstraints()
{
	const float DesiredHorizontalDistance = HorizontalDistance;
	const float DesiredVerticalDistance = VerticalDistance;

	for (int32 Iteration = 0; Iteration < SolverIterationCount; ++Iteration)
	{
		// Solve horizontal constraint
		for (int32 Y = 0; Y < VerticalVertexCount; ++Y)
		{
			for (int32 X = 0; X < HorizontalVertexCount - 1; ++X)
			{
				FClothParticle& ParticleA = ClothParticles[Y * HorizontalVertexCount + X];
				FClothParticle& ParticleB = ClothParticles[Y * HorizontalVertexCount + X + 1];

				SolveDistanceConstraint(ParticleA, ParticleB, DesiredHorizontalDistance);

#if VPS_CLOTH_DEBUG
				SolveSphericalConstraint(ParticleA, SphereRadius, SphereLoc);
				SolveSphericalConstraint(ParticleA, SphereRadius, SphereLoc);
				SolveGroundConstraint(ParticleA, GroundZ);
				SolveGroundConstraint(ParticleB, GroundZ);
#endif
			}
		}

		// Solve vertical constraint
		for (int32 X = 0; X < HorizontalVertexCount; ++X)
		{
			for (int32 Y = 0; Y < VerticalVertexCount - 1; ++Y)
			{
				FClothParticle& ParticleA = ClothParticles[Y * HorizontalVertexCount + X];
				FClothParticle& ParticleB = ClothParticles[(Y + 1) * HorizontalVertexCount + X];

				SolveDistanceConstraint(ParticleA, ParticleB, DesiredVerticalDistance);

#if VPS_CLOTH_DEBUG
				SolveSphericalConstraint(ParticleA, SphereRadius, SphereLoc);
				SolveSphericalConstraint(ParticleA, SphereRadius, SphereLoc);
				SolveGroundConstraint(ParticleA, GroundZ);
				SolveGroundConstraint(ParticleB, GroundZ);
#endif
			}
		}
	}

	// Copy to position array
	ParallelFor(ClothParticles.Num(), [&](int32 Index)
	{
		ClothPositions[Index] = ClothParticles[Index].NewLoc;
	});
}

void UProceduralClothComponent::ComputeNormals()
{
	ParallelFor(ClothParticles.Num(), [&](int32 Index)
	{
		const int32 X = Index % HorizontalVertexCount;
		const int32 Y = Index / HorizontalVertexCount;
		const FVector ParticlePos = ClothPositions[Index];
		FVector DirA, DirB;
		int32 IndexA, IndexB;

		IndexA = X > 0 ? X - 1 : 1; // Favor left side
		IndexB = Y > 0 ? Y - 1 : 1; // Favor up side

		FVector PosA = ClothPositions[Y * HorizontalVertexCount + IndexA]; // Left side position
		FVector PosB = ClothPositions[IndexB * HorizontalVertexCount + X]; // Up side position

		DirA = PosA - ParticlePos;
		DirB = PosB - ParticlePos;

		ClothNormals[Index] = (DirB ^ DirA).GetSafeNormal();

		if ((X == 0) ^ (Y == 0))
		{
			ClothNormals[Index] *= -1;
		}
	});
}

void UProceduralClothComponent::PerformSubstepParallel(float InSubstepTime, const FVector& Gravity)
{
	VerletIntegrateParallel(InSubstepTime, Gravity);
	SolveConstraintsParallel();
	ComputeNormalsParallel();
}

void UProceduralClothComponent::VerletIntegrateParallel(float InSubstepTime, const FVector& Gravity)
{
	ENQUEUE_RENDER_COMMAND(ClothComputeCommands)
	(
		[this, InSubstepTime, Gravity](FRHICommandListImmediate& RHICmdList)
		{
			check(IsInRenderingThread());

			TShaderMapRef<FClothVerletComputeShader> ClothVerletComputeShader(GetGlobalShaderMap(ERHIFeatureLevel::SM5));
			RHICmdList.SetComputeShader(ClothVerletComputeShader.GetComputeShader());

			// Bind shader buffers
			ClothVerletComputeShader->BindShaderBuffers(RHICmdList, ParticlesStructuredBufferUAV);

			// Bind shader uniform
			FClothVerletComputeShaderParameters UniformParam;
			UniformParam.SizeX = HorizontalVertexCount;
			UniformParam.SizeY = VerticalVertexCount;
			UniformParam.SubStepSqr = InSubstepTime * InSubstepTime;
			UniformParam.Gravity = Gravity;
			ClothVerletComputeShader->SetShaderParameters(RHICmdList, UniformParam);

			// Dispatch shader
			const int ThreadGroupCountX = FMath::CeilToInt(HorizontalVertexCount / 32.f);
			const int ThreadGroupCountY = FMath::CeilToInt(VerticalVertexCount / 32.f);
			const int ThreadGroupCountZ = 1;
			DispatchComputeShader(RHICmdList, ClothVerletComputeShader, ThreadGroupCountX, ThreadGroupCountY, ThreadGroupCountZ);

			// Unbind shader buffers
			ClothVerletComputeShader->UnbindShaderBuffers(RHICmdList);
		}
	);
}

void UProceduralClothComponent::SolveConstraintsParallel()
{
	ENQUEUE_RENDER_COMMAND(ClothComputeCommands)
	(
		[this](FRHICommandListImmediate& RHICmdList)
		{
			check(IsInRenderingThread());

			TShaderMapRef<FClothConstraintComputeShader> ClothConstraintComputeShader(GetGlobalShaderMap(ERHIFeatureLevel::SM5));
			RHICmdList.SetComputeShader(ClothConstraintComputeShader.GetComputeShader());

			for (int32 Iteration = 0; Iteration < SolverIterationCount; ++Iteration)
			{
				for (uint32 Direction = 0; Direction <= 1; ++Direction)
				{
					// Make UAV safe for read
					RHICmdList.TransitionResource(
						EResourceTransitionAccess::ERWBarrier,
						EResourceTransitionPipeline::EComputeToCompute,
						PositionsStructuredBufferUAV,
						nullptr);

					// Bind shader buffers
					ClothConstraintComputeShader->BindShaderBuffers(RHICmdList, ParticlesStructuredBufferUAV, PositionsStructuredBufferUAV);

					// Bind shader uniform
					FClothConstraintComputeShaderParameters UniformParam;
					UniformParam.SizeX           = HorizontalVertexCount;
					UniformParam.SizeY           = VerticalVertexCount;
					UniformParam.Direction       = Direction;
					UniformParam.GroundZ         = GroundZ;
					UniformParam.SphereRadius    = SphereRadius;
					UniformParam.SphereLoc       = SphereLoc;
					UniformParam.DesiredDistance = (Direction == 0) ? HorizontalDistance : VerticalDistance;
					ClothConstraintComputeShader->SetShaderParameters(RHICmdList, UniformParam);

					// Dispatch shader
					const int ThreadGroupCountX = FMath::CeilToInt(((Direction == 0) ? VerticalVertexCount : HorizontalVertexCount) / 64.f);
					const int ThreadGroupCountY = 1;
					const int ThreadGroupCountZ = 1;
					DispatchComputeShader(RHICmdList, ClothConstraintComputeShader, ThreadGroupCountX, ThreadGroupCountY, ThreadGroupCountZ);

					// Unbind shader buffers
					ClothConstraintComputeShader->UnbindShaderBuffers(RHICmdList);
				}
			}

			// Make UAV safe for read
			RHICmdList.TransitionResource(
				EResourceTransitionAccess::ERWBarrier,
				EResourceTransitionPipeline::EComputeToCompute,
				PositionsStructuredBufferUAV,
				nullptr);

			// Read back the position buffer. No need to copy back the particle buffer since all we care about is position
			const int32 BufferSize = HorizontalVertexCount * VerticalVertexCount * sizeof(FVector);
			FVector* SrcPtr = (FVector*)RHILockStructuredBuffer(PositionsStructuredBuffer.GetReference(), 0, BufferSize, EResourceLockMode::RLM_ReadOnly);
			FVector* DstPtr = ClothPositions.GetData();
			FMemory::Memcpy(DstPtr, SrcPtr, BufferSize);
			RHIUnlockStructuredBuffer(PositionsStructuredBuffer.GetReference());
		}
	);
}

void UProceduralClothComponent::ComputeNormalsParallel()
{
	ENQUEUE_RENDER_COMMAND(ClothComputeCommands)
	(
		[this](FRHICommandListImmediate& RHICmdList)
		{
			check(IsInRenderingThread());

			// Make UAV safe for read
			RHICmdList.TransitionResource(
				EResourceTransitionAccess::ERWBarrier,
				EResourceTransitionPipeline::EComputeToCompute,
				PositionsStructuredBufferUAV,
				nullptr);

			TShaderMapRef<FClothNormalComputeShader> ClothNormalComputeShader(GetGlobalShaderMap(ERHIFeatureLevel::SM5));
			RHICmdList.SetComputeShader(ClothNormalComputeShader.GetComputeShader());

			// Bind shader buffers
			ClothNormalComputeShader->BindShaderBuffers(RHICmdList, PositionsStructuredBufferUAV, NormalsStructuredBufferUAV);

			// Bind shader uniform
			FClothNormalComputeShaderParameters UniformParam;
			UniformParam.SizeX = HorizontalVertexCount;
			UniformParam.SizeY = VerticalVertexCount;
			ClothNormalComputeShader->SetShaderParameters(RHICmdList, UniformParam);

			// Dispatch shader
			const int ThreadGroupCountX = FMath::CeilToInt(HorizontalVertexCount / 32.f);
			const int ThreadGroupCountY = FMath::CeilToInt(VerticalVertexCount / 32.f);
			const int ThreadGroupCountZ = 1;
			DispatchComputeShader(RHICmdList, ClothNormalComputeShader, ThreadGroupCountX, ThreadGroupCountY, ThreadGroupCountZ);

			// Unbind shader buffers
			ClothNormalComputeShader->UnbindShaderBuffers(RHICmdList);

			// Make UAV safe for read
			RHICmdList.TransitionResource(
				EResourceTransitionAccess::ERWBarrier,
				EResourceTransitionPipeline::EComputeToCompute,
				NormalsStructuredBufferUAV,
				nullptr);

			// Read back the normal buffer
			const int32 BufferSize = HorizontalVertexCount * VerticalVertexCount * sizeof(FVector);
			FVector* SrcPtr = (FVector*)RHILockStructuredBuffer(NormalsStructuredBuffer.GetReference(), 0, BufferSize, EResourceLockMode::RLM_ReadOnly);
			FVector* DstPtr = ClothNormals.GetData();
			FMemory::Memcpy(DstPtr, SrcPtr, BufferSize);
			RHIUnlockStructuredBuffer(NormalsStructuredBuffer.GetReference());
		}
	);
}

void UProceduralClothComponent::ReleaseBufferResources()
{
	ReleaseBufferResource<FStructuredBufferRHIRef>(ParticlesStructuredBuffer);
	ReleaseBufferResource<FStructuredBufferRHIRef>(PositionsStructuredBuffer);
	ReleaseBufferResource<FStructuredBufferRHIRef>(NormalsStructuredBuffer);
	ReleaseBufferResource<FUnorderedAccessViewRHIRef>(ParticlesStructuredBufferUAV);
	ReleaseBufferResource<FUnorderedAccessViewRHIRef>(PositionsStructuredBufferUAV);
	ReleaseBufferResource<FUnorderedAccessViewRHIRef>(NormalsStructuredBufferUAV);
}

void UProceduralClothComponent::TickComponent(float DeltaTime, enum ELevelTick TickType, FActorComponentTickFunction *ThisTickFunction)
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

	// Tick cloth step and solver
	TickClothComponent(DeltaTime);

	// Update the mesh vertex positions after solver
	UpdateProceduralMesh(true);

	// Need to send new data to render thread
	if (CastShadow)
	{
		MarkRenderStateDirty();
	}
	else
	{
		MarkRenderDynamicDataDirty();
	}

	// Call this because bounds have changed
	UpdateComponentToWorld();
}

void UProceduralClothComponent::TickClothComponent(float DeltaTime)
{
	const FVector Gravity = FVector(0, 0, GetWorld()->GetGravityZ()) * ClothGravityScale;

	// Ensure a non-zero substep
	float UseSubstep = FMath::Max(SubstepTime, 0.005f);

	// Perform simulation substeps
	TimeRemainder += DeltaTime;
	while (TimeRemainder > UseSubstep)
	{
		if (bUseGPUAcceleration)
		{
			PerformSubstepParallel(UseSubstep, Gravity);
		}
		else
		{
			PerformSubstep(UseSubstep, Gravity);
		}
		TimeRemainder -= UseSubstep;
	}
}
