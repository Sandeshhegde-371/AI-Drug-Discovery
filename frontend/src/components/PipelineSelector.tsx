import { Sparkles, Search, BarChart3, ListOrdered, Check } from 'lucide-react';
import { Card } from './ui/card';
import { Badge } from './ui/badge';

export type Pipeline = 'generative' | 'virtual' | 'properties' | 'prioritize';

interface PipelineSelectorProps {
  selectedPipeline: Pipeline | null;
  onSelectPipeline: (pipeline: Pipeline) => void;
}

const pipelines = [
  {
    id: 'generative' as Pipeline,
    name: 'Generative Screening',
    icon: Sparkles,
    description: 'Generate novel molecules based on your seed data using AI molecular design',
    features: ['Create new candidates', 'Optimize drug-like properties', 'Explore chemical space'],
    gradient: 'from-purple-600 to-indigo-600',
    bgGradient: 'from-purple-50 to-indigo-50',
    borderColor: 'border-purple-300',
    hoverBorder: 'hover:border-purple-400',
  },
  {
    id: 'virtual' as Pipeline,
    name: 'Virtual Screening',
    icon: Search,
    description: 'Screen large molecular libraries to identify the most promising candidates',
    features: ['Fast compound filtering', 'Target binding prediction', 'Similarity scoring'],
    gradient: 'from-blue-600 to-cyan-600',
    bgGradient: 'from-blue-50 to-cyan-50',
    borderColor: 'border-blue-300',
    hoverBorder: 'hover:border-blue-400',
  },
  {
    id: 'properties' as Pipeline,
    name: 'Predict Properties',
    icon: BarChart3,
    description: 'Predict physicochemical and biological properties for your molecules',
    features: ['Activity prediction', 'ADMET properties', 'Toxicity assessment'],
    gradient: 'from-green-600 to-emerald-600',
    bgGradient: 'from-green-50 to-emerald-50',
    borderColor: 'border-green-300',
    hoverBorder: 'hover:border-green-400',
  },
  {
    id: 'prioritize' as Pipeline,
    name: 'Prioritize Candidates',
    icon: ListOrdered,
    description: 'Rank molecules based on multiple criteria to optimize lab testing order',
    features: ['Multi-property scoring', 'Custom ranking criteria', 'Resource optimization'],
    gradient: 'from-orange-600 to-amber-600',
    bgGradient: 'from-orange-50 to-amber-50',
    borderColor: 'border-orange-300',
    hoverBorder: 'hover:border-orange-400',
  },
];

export function PipelineSelector({ selectedPipeline, onSelectPipeline }: PipelineSelectorProps) {
  return (
    <Card className="p-8 bg-white/70 backdrop-blur-sm border-purple-100 shadow-lg shadow-purple-100/50">
      <div className="mb-8">
        <div className="flex items-center gap-2 mb-3">
          <div className="w-10 h-10 bg-gradient-to-br from-purple-600 to-indigo-600 rounded-lg flex items-center justify-center">
            <Sparkles className="w-5 h-5 text-white" />
          </div>
          <div>
            <h3 className="text-gray-900">Select AI Pipeline</h3>
            <p className="text-sm text-gray-600">
              Choose the analysis that fits your research goals
            </p>
          </div>
        </div>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-2 gap-5">
        {pipelines.map((pipeline) => {
          const Icon = pipeline.icon;
          const isSelected = selectedPipeline === pipeline.id;
          
          return (
            <div
              key={pipeline.id}
              onClick={() => onSelectPipeline(pipeline.id)}
              className={`relative group border-2 rounded-xl p-6 cursor-pointer transition-all duration-300 ${
                isSelected
                  ? `${pipeline.borderColor} bg-gradient-to-br ${pipeline.bgGradient} shadow-lg scale-[1.02]`
                  : `border-gray-200 bg-white ${pipeline.hoverBorder} hover:shadow-md hover:scale-[1.01]`
              }`}
            >
              {/* Selection indicator */}
              {isSelected && (
                <div className="absolute top-4 right-4">
                  <div className={`w-8 h-8 rounded-full bg-gradient-to-br ${pipeline.gradient} flex items-center justify-center`}>
                    <Check className="w-5 h-5 text-white" />
                  </div>
                </div>
              )}

              <div className="flex items-start gap-4 mb-4">
                <div className={`w-12 h-12 rounded-xl bg-gradient-to-br ${pipeline.gradient} flex items-center justify-center flex-shrink-0 shadow-md`}>
                  <Icon className="w-6 h-6 text-white" />
                </div>
                <div className="flex-1">
                  <h4 className="text-gray-900 mb-1">{pipeline.name}</h4>
                  <p className="text-sm text-gray-600 leading-relaxed">{pipeline.description}</p>
                </div>
              </div>
              
              <div className="space-y-2 pl-16">
                {pipeline.features.map((feature, index) => (
                  <div key={index} className="flex items-center gap-2 text-sm text-gray-700">
                    <div className={`w-1.5 h-1.5 rounded-full ${isSelected ? 'bg-purple-600' : 'bg-gray-400'}`} />
                    <span>{feature}</span>
                  </div>
                ))}
              </div>
            </div>
          );
        })}
      </div>
    </Card>
  );
}
