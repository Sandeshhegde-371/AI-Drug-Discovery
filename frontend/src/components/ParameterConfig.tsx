import { useState } from 'react';
import { Settings, Play, Sparkles } from 'lucide-react';
import { Card } from './ui/card';
import { Button } from './ui/button';
import { Label } from './ui/label';
import { Input } from './ui/input';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from './ui/select';
import { Slider } from './ui/slider';
import type { Pipeline } from './PipelineSelector';

interface ParameterConfigProps {
  pipeline: Pipeline;
  onRunPipeline: (params: any) => void;
  isRunning: boolean;
}

export function ParameterConfig({ pipeline, onRunPipeline, isRunning }: ParameterConfigProps) {
  const [params, setParams] = useState<any>({
    numCandidates: 100,
    similarityThreshold: 0.7,
    propertyFilter: 'all',
    rankingCriteria: 'activity',
    confidenceThreshold: 0.8,
  });

  const updateParam = (key: string, value: any) => {
    setParams({ ...params, [key]: value });
  };

  const renderGenerativeParams = () => (
    <>
      <div className="space-y-3">
        <Label htmlFor="numCandidates" className="text-gray-900">Number of Candidates to Generate</Label>
        <Input
          id="numCandidates"
          type="number"
          value={params.numCandidates}
          onChange={(e) => updateParam('numCandidates', parseInt(e.target.value))}
          min="10"
          max="1000"
          className="border-purple-200 focus:border-purple-400"
        />
        <p className="text-sm text-gray-500">Generate between 10 and 1000 novel molecules</p>
      </div>
      
      <div className="space-y-3">
        <div className="flex items-center justify-between">
          <Label className="text-gray-900">Similarity to Seed Molecules</Label>
          <span className="text-sm font-semibold text-purple-600">{(params.similarityThreshold * 100).toFixed(0)}%</span>
        </div>
        <Slider
          value={[params.similarityThreshold * 100]}
          onValueChange={(val) => updateParam('similarityThreshold', val[0] / 100)}
          min={30}
          max={95}
          step={5}
          className="[&_[role=slider]]:bg-purple-600 [&_[role=slider]]:border-purple-600"
        />
        <p className="text-sm text-gray-500">Higher values generate molecules more similar to your input</p>
      </div>

      <div className="space-y-3">
        <Label htmlFor="drugLikeness" className="text-gray-900">Drug-Likeness Filter</Label>
        <Select value={params.propertyFilter} onValueChange={(val) => updateParam('propertyFilter', val)}>
          <SelectTrigger id="drugLikeness" className="border-purple-200">
            <SelectValue />
          </SelectTrigger>
          <SelectContent>
            <SelectItem value="all">All molecules</SelectItem>
            <SelectItem value="lipinski">Lipinski's Rule of Five</SelectItem>
            <SelectItem value="veber">Veber Rules</SelectItem>
            <SelectItem value="custom">Custom filters</SelectItem>
          </SelectContent>
        </Select>
      </div>
    </>
  );

  const renderVirtualParams = () => (
    <>
      <div className="space-y-3">
        <Label htmlFor="screeningSize" className="text-gray-900">Library Size to Screen</Label>
        <Select value={params.numCandidates.toString()} onValueChange={(val) => updateParam('numCandidates', parseInt(val))}>
          <SelectTrigger id="screeningSize" className="border-blue-200">
            <SelectValue />
          </SelectTrigger>
          <SelectContent>
            <SelectItem value="1000">1,000 compounds</SelectItem>
            <SelectItem value="10000">10,000 compounds</SelectItem>
            <SelectItem value="100000">100,000 compounds</SelectItem>
            <SelectItem value="1000000">1,000,000 compounds</SelectItem>
          </SelectContent>
        </Select>
      </div>

      <div className="space-y-3">
        <div className="flex items-center justify-between">
          <Label className="text-gray-900">Binding Affinity Threshold</Label>
          <span className="text-sm font-semibold text-blue-600">{(params.confidenceThreshold * 100).toFixed(0)}%</span>
        </div>
        <Slider
          value={[params.confidenceThreshold * 100]}
          onValueChange={(val) => updateParam('confidenceThreshold', val[0] / 100)}
          min={50}
          max={95}
          step={5}
          className="[&_[role=slider]]:bg-blue-600 [&_[role=slider]]:border-blue-600"
        />
        <p className="text-sm text-gray-500">Minimum predicted binding affinity to report</p>
      </div>

      <div className="space-y-3">
        <Label htmlFor="scoringMethod" className="text-gray-900">Scoring Method</Label>
        <Select value={params.rankingCriteria} onValueChange={(val) => updateParam('rankingCriteria', val)}>
          <SelectTrigger id="scoringMethod" className="border-blue-200">
            <SelectValue />
          </SelectTrigger>
          <SelectContent>
            <SelectItem value="activity">Predicted Activity</SelectItem>
            <SelectItem value="docking">Docking Score</SelectItem>
            <SelectItem value="similarity">Structural Similarity</SelectItem>
            <SelectItem value="composite">Composite Score</SelectItem>
          </SelectContent>
        </Select>
      </div>
    </>
  );

  const renderPropertiesParams = () => (
    <>
      <div className="space-y-3">
        <Label className="text-gray-900">Properties to Predict</Label>
        <div className="space-y-2 p-5 border-2 border-green-200 bg-green-50/30 rounded-lg">
          {['Activity', 'Solubility', 'Toxicity', 'Permeability', 'Metabolic Stability', 'hERG Inhibition'].map((prop) => (
            <label key={prop} className="flex items-center gap-3 cursor-pointer group">
              <input type="checkbox" defaultChecked className="rounded border-green-300 text-green-600 focus:ring-green-500" />
              <span className="text-gray-700 group-hover:text-gray-900">{prop}</span>
            </label>
          ))}
        </div>
      </div>

      <div className="space-y-3">
        <div className="flex items-center justify-between">
          <Label className="text-gray-900">Prediction Confidence</Label>
          <span className="text-sm font-semibold text-green-600">{(params.confidenceThreshold * 100).toFixed(0)}%</span>
        </div>
        <Slider
          value={[params.confidenceThreshold * 100]}
          onValueChange={(val) => updateParam('confidenceThreshold', val[0] / 100)}
          min={60}
          max={95}
          step={5}
          className="[&_[role=slider]]:bg-green-600 [&_[role=slider]]:border-green-600"
        />
        <p className="text-sm text-gray-500">Minimum confidence level for predictions</p>
      </div>
    </>
  );

  const renderPrioritizeParams = () => (
    <>
      <div className="space-y-3">
        <Label htmlFor="rankingCriteria" className="text-gray-900">Primary Ranking Criterion</Label>
        <Select value={params.rankingCriteria} onValueChange={(val) => updateParam('rankingCriteria', val)}>
          <SelectTrigger id="rankingCriteria" className="border-orange-200">
            <SelectValue />
          </SelectTrigger>
          <SelectContent>
            <SelectItem value="activity">Predicted Activity</SelectItem>
            <SelectItem value="safety">Safety Profile</SelectItem>
            <SelectItem value="druglikeness">Drug-Likeness Score</SelectItem>
            <SelectItem value="synthesizability">Ease of Synthesis</SelectItem>
            <SelectItem value="balanced">Balanced Multi-Criteria</SelectItem>
          </SelectContent>
        </Select>
      </div>

      <div className="space-y-3">
        <Label htmlFor="topN" className="text-gray-900">Number of Top Candidates</Label>
        <Input
          id="topN"
          type="number"
          value={params.numCandidates}
          onChange={(e) => updateParam('numCandidates', parseInt(e.target.value))}
          min="10"
          max="500"
          className="border-orange-200 focus:border-orange-400"
        />
        <p className="text-sm text-gray-500">How many prioritized candidates to return</p>
      </div>

      <div className="space-y-3">
        <div className="flex items-center justify-between">
          <Label className="text-gray-900">Diversity Constraint</Label>
          <span className="text-sm font-semibold text-orange-600">{(params.similarityThreshold * 100).toFixed(0)}%</span>
        </div>
        <Slider
          value={[params.similarityThreshold * 100]}
          onValueChange={(val) => updateParam('similarityThreshold', val[0] / 100)}
          min={30}
          max={90}
          step={5}
          className="[&_[role=slider]]:bg-orange-600 [&_[role=slider]]:border-orange-600"
        />
        <p className="text-sm text-gray-500">Higher values ensure more structural diversity in results</p>
      </div>
    </>
  );

  const renderParams = () => {
    switch (pipeline) {
      case 'generative':
        return renderGenerativeParams();
      case 'virtual':
        return renderVirtualParams();
      case 'properties':
        return renderPropertiesParams();
      case 'prioritize':
        return renderPrioritizeParams();
      default:
        return null;
    }
  };

  return (
    <Card className="p-8 bg-white/70 backdrop-blur-sm border-purple-100 shadow-lg shadow-purple-100/50">
      <div className="mb-8">
        <div className="flex items-center gap-2 mb-3">
          <div className="w-10 h-10 bg-gradient-to-br from-purple-600 to-indigo-600 rounded-lg flex items-center justify-center">
            <Settings className="w-5 h-5 text-white" />
          </div>
          <div>
            <h3 className="text-gray-900">Configure Parameters</h3>
            <p className="text-sm text-gray-600">
              Customize settings for your analysis
            </p>
          </div>
        </div>
      </div>

      <div className="space-y-6">
        {renderParams()}
      </div>

      <div className="mt-8 pt-6 border-t border-gray-200">
        <Button
          onClick={() => onRunPipeline(params)}
          disabled={isRunning}
          className="w-full gap-2 h-12 bg-gradient-to-r from-purple-600 to-indigo-600 hover:from-purple-700 hover:to-indigo-700 shadow-lg shadow-purple-500/30"
          size="lg"
        >
          {isRunning ? (
            <>
              <div className="w-5 h-5 border-2 border-white border-t-transparent rounded-full animate-spin" />
              <span>Processing Pipeline...</span>
            </>
          ) : (
            <>
              <Play className="w-5 h-5" />
              <span>Run AI Pipeline</span>
              <Sparkles className="w-4 h-4 ml-1" />
            </>
          )}
        </Button>
      </div>
    </Card>
  );
}
