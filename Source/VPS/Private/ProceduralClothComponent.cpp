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
	SHADER_PARAMETER(int,   SizeX)
	SHADER_PARAMETER(int,   SizeY)
	SHADER_PARAMETER(float, DesiredDistance)
	SHADER_PARAMETER(int,   Direction)   // 0 for solving horizontal constraint. 1 for solving vertical constraint
END_GLOBAL_SHADER_PARAMETER_STRUCT()
IMPLEMENT_GLOBAL_SHADER_PARAMETER_STRUCT(FClothConstraintComputeShaderParameters, "ClothConstraintUniform");

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

	virtual bool Serialize(FArchive& Ar) override
	{
		bool bShaderHasOutdatedParameters = FGlobalShader::Serialize(Ar);
		Ar << OutputParticleBuffer;
		return bShaderHasOutdatedParameters;
	}

	void BindShaderBuffers(FRHICommandList& RHICmdList, FUnorderedAccessViewRHIRef OutputBufferUAV)
	{
		FRHIComputeShader* ComputeShaderRHI = GetComputeShader();
		SetUAVParameter(RHICmdList, ComputeShaderRHI, OutputParticleBuffer, OutputBufferUAV);
	}

	void UnbindShaderBuffers(FRHICommandList& RHICmdList)
	{
		FRHIComputeShader* ComputeShaderRHI = GetComputeShader();
		SetUAVParameter(RHICmdList, ComputeShaderRHI, OutputParticleBuffer, FUnorderedAccessViewRHIRef());
	}

	void SetShaderParameters(FRHICommandList& RHICmdList, const FClothVerletComputeShaderParameters& Parameters)
	{
		FRHIComputeShader* ComputeShaderRHI = GetComputeShader();
		SetUniformBufferParameterImmediate(RHICmdList, ComputeShaderRHI, GetUniformBufferParameter<FClothVerletComputeShaderParameters>(), Parameters);
	}

private:

	FShaderResourceParameter OutputParticleBuffer;
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
	}

	static bool ShouldCompilePermutation(const FGlobalShaderPermutationParameters& Parameters)
	{
		return IsFeatureLevelSupported(Parameters.Platform, ERHIFeatureLevel::SM5);
	}

	static bool ShouldCache(EShaderPlatform Platform)
	{
		return IsFeatureLevelSupported(Platform, ERHIFeatureLevel::SM5);
	}

	virtual bool Serialize(FArchive& Ar) override
	{
		bool bShaderHasOutdatedParameters = FGlobalShader::Serialize(Ar);
		Ar << OutputParticleBuffer;
		return bShaderHasOutdatedParameters;
	}

	void BindShaderBuffers(FRHICommandList& RHICmdList, FUnorderedAccessViewRHIRef OutputBufferUAV)
	{
		FRHIComputeShader* ComputeShaderRHI = GetComputeShader();
		SetUAVParameter(RHICmdList, ComputeShaderRHI, OutputParticleBuffer, OutputBufferUAV);
	}

	void UnbindShaderBuffers(FRHICommandList& RHICmdList)
	{
		FRHIComputeShader* ComputeShaderRHI = GetComputeShader();
		SetUAVParameter(RHICmdList, ComputeShaderRHI, OutputParticleBuffer, FUnorderedAccessViewRHIRef());
	}

	void SetShaderParameters(FRHICommandList& RHICmdList, const FClothConstraintComputeShaderParameters& Parameters)
	{
		FRHIComputeShader* ComputeShaderRHI = GetComputeShader();
		SetUniformBufferParameterImmediate(RHICmdList, ComputeShaderRHI, GetUniformBufferParameter<FClothConstraintComputeShaderParameters>(), Parameters);
	}

private:

	FShaderResourceParameter OutputParticleBuffer;
};

IMPLEMENT_SHADER_TYPE(, FClothVerletComputeShader, TEXT("/Plugin/VPS/ClothConstraintComputeShader.usf"), TEXT("ComputeClothVerlet"), SF_Compute);
IMPLEMENT_SHADER_TYPE(, FClothConstraintComputeShader, TEXT("/Plugin/VPS/ClothConstraintComputeShader.usf"), TEXT("SolveClothConstraint"), SF_Compute);

namespace
{
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

	PrimaryComponentTick.bCanEverTick = true;
	PrimaryComponentTick.bStartWithTickEnabled = true;
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

		FRHIResourceCreateInfo CreateInfo;
		CreateInfo.ResourceArray = &ClothParticles;

		ParticlesStructuredBuffer = RHICreateStructuredBuffer(
			sizeof(FClothParticle),                   // Stride
			sizeof(FClothParticle) * VertexCount,     // Size
			BUF_UnorderedAccess | BUF_ShaderResource, // Usage
			CreateInfo                                // Create info
		);

		ParticlesStructuredBufferUAV = RHICreateUnorderedAccessView(ParticlesStructuredBuffer, true, false);
	}
}

void UProceduralClothComponent::InitClothParticles()
{
	const int32 VertexCount = HorizontalVertexCount * VerticalVertexCount;

	ClothParticles.Reset();
	ClothParticles.AddUninitialized(VertexCount);

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
			ClothParticle.bFree = (X - HorizontalVertexCount + 1) != 0 && (Y - VerticalVertexCount + 1) != 0;
		}
	}
}

void UProceduralClothComponent::UpdateProceduralMesh(bool bVertexOnly)
{
	const int32 VertexCount = HorizontalVertexCount * VerticalVertexCount;

	TArray<FVector> Vertices;
	TArray<FVector> Normals;
	TArray<FVector2D> UVs;
	TArray<FColor> VertexColors;
	TArray<FProcMeshTangent> Tangents;

	Vertices.AddUninitialized(VertexCount);

	ParallelFor(Vertices.Num(), [&](int32 Index)
	{
		Vertices[Index] = ClothParticles[Index].NewLoc;
	});

	if (!bVertexOnly)
	{
		const int32 TriangleCount = (HorizontalVertexCount - 1) * (VerticalVertexCount - 1) * 6;
		int32 Idx = 0;
		TArray<int32> Triangles;

		Normals.AddUninitialized(VertexCount);
		UVs.AddUninitialized(VertexCount);
		Triangles.AddUninitialized(TriangleCount);

		for (int32 Y = 0; Y < VerticalVertexCount; ++Y)
		{
			for (int32 X = 0; X < HorizontalVertexCount; ++X)
			{
				float U = StaticCast<float>(X) / HorizontalVertexCount - 1;
				float V = StaticCast<float>(Y) / VerticalVertexCount - 1;
				int32 Index = Y * HorizontalVertexCount + X;

				Normals[Index] = FVector::UpVector;
				UVs[Index] = FVector2D(U, V);
			}
		}

		for (int32 Y = 0; Y < VerticalVertexCount - 1; ++Y)
		{
			for (int32 X = 0; X < HorizontalVertexCount - 1; ++X)
			{
				int A = Y * HorizontalVertexCount + X;
				int B = A + HorizontalVertexCount;
				int C = A + HorizontalVertexCount + 1;
				int D = A + 1;

				Triangles[Idx++] = A;
				Triangles[Idx++] = B;
				Triangles[Idx++] = C;

				Triangles[Idx++] = A;
				Triangles[Idx++] = C;
				Triangles[Idx++] = D;
			}
		}

		ClearMeshSection(0);
		CreateMeshSection(0, Vertices, Triangles, Normals, UVs, VertexColors, Tangents, false);
	}

	else
	{
		UpdateMeshSection(0, Vertices, Normals, UVs, VertexColors, Tangents);
	}
}

void UProceduralClothComponent::PerformSubstep(float InSubstepTime, const FVector& Gravity)
{
	VerletIntegrate(InSubstepTime, Gravity);
	SolveConstraints();
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
			}
		}
	}
}

void UProceduralClothComponent::PerformSubstepParallel(float InSubstepTime, const FVector& Gravity)
{
	VerletIntegrateParallel(InSubstepTime, Gravity);
	SolveConstraintsParallel();
}

void UProceduralClothComponent::VerletIntegrateParallel(float InSubstepTime, const FVector& Gravity)
{
	ENQUEUE_RENDER_COMMAND(ClothComputeCommands)
	(
		[this, InSubstepTime, Gravity](FRHICommandListImmediate& RHICmdList)
		{
			check(IsInRenderingThread());

			TShaderMapRef<FClothVerletComputeShader> ClothVerletComputeShader(GetGlobalShaderMap(ERHIFeatureLevel::SM5));
			RHICmdList.SetComputeShader(ClothVerletComputeShader->GetComputeShader());

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
			const int ThreadGroupCountX = StaticCast<int>(HorizontalVertexCount / 32);
			const int ThreadGroupCountY = StaticCast<int>(VerticalVertexCount / 32);
			const int ThreadGroupCountZ = 1;
			DispatchComputeShader(RHICmdList, *ClothVerletComputeShader, ThreadGroupCountX, ThreadGroupCountY, ThreadGroupCountZ);

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
			RHICmdList.SetComputeShader(ClothConstraintComputeShader->GetComputeShader());

			for (int32 Iteration = 0; Iteration < SolverIterationCount; ++Iteration)
			{
				for (uint32 Direction = 0; Direction <= 1; ++Direction)
				{
					// Bind shader buffers
					ClothConstraintComputeShader->BindShaderBuffers(RHICmdList, ParticlesStructuredBufferUAV);

					// Bind shader uniform
					FClothConstraintComputeShaderParameters UniformParam;
					UniformParam.SizeX = HorizontalVertexCount;
					UniformParam.SizeY = VerticalVertexCount;
					UniformParam.Direction = Direction;
					UniformParam.DesiredDistance = (Direction == 0) ? HorizontalDistance : VerticalDistance;
					ClothConstraintComputeShader->SetShaderParameters(RHICmdList, UniformParam);

					// Dispatch shader
					const int ThreadGroupCountX = StaticCast<int>(((Direction == 0) ? VerticalVertexCount : HorizontalVertexCount) / 64);
					const int ThreadGroupCountY = 1;
					const int ThreadGroupCountZ = 1;
					DispatchComputeShader(RHICmdList, *ClothConstraintComputeShader, ThreadGroupCountX, ThreadGroupCountY, ThreadGroupCountZ);

					// Unbind shader buffers
					ClothConstraintComputeShader->UnbindShaderBuffers(RHICmdList);
				}
			}

			// Read back
			const int32 BufferSize = HorizontalVertexCount * VerticalVertexCount * sizeof(FClothParticle);
			FClothParticle* SrcPtr = (FClothParticle*) RHILockStructuredBuffer(ParticlesStructuredBuffer.GetReference(), 0, BufferSize, EResourceLockMode::RLM_ReadOnly);
			FClothParticle* DstPtr = ClothParticles.GetData();
			FMemory::Memcpy(DstPtr, SrcPtr, BufferSize);
			RHIUnlockStructuredBuffer(ParticlesStructuredBuffer.GetReference());
		}
	);
}

void UProceduralClothComponent::ReleaseBufferResources()
{
#define SafeReleaseBufferResource(Texture)   \
	do {                                     \
		if (Texture.IsValid()) {             \
			Texture->Release();              \
		}                                    \
	} while(0);

	SafeReleaseBufferResource(ParticlesStructuredBuffer);
	SafeReleaseBufferResource(ParticlesStructuredBufferUAV);
}

void UProceduralClothComponent::TickComponent(float DeltaTime, enum ELevelTick TickType, FActorComponentTickFunction *ThisTickFunction)
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

	// Tick cloth step and solver
	TickClothComponent(DeltaTime);

	// Update the mesh vertex positions after solver
	UpdateProceduralMesh(true);

	// Need to send new data to render thread
	MarkRenderDynamicDataDirty();

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